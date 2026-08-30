//
// Harmonisation stage 2: merge multiple callers' already-clustered site
// VCFs (e.g. Manta's and Wham's output from vcfs_cluster_svcluster) into
// one cross-caller site list. Mirrors GATK-SV's CombineBatches workflow
// (wdl/CombineBatches.wdl) -- see
// README.md#how-gatk-svs-harmonisation-actually-works for the stage-1
// vs stage-2 distinction.
//
// Steps, matching CombineBatches:
//   1. Naive join across callers (SVCluster, near-zero overlap
//      thresholds -- only merges exact-position duplicates)
//   2. Real clustering pass (SVCluster, realistic overlap thresholds)
//   3+4. Two rounds of context-aware re-clustering (GroupedSVCluster,
//      stratified by genomic context tracks: simple repeats, segmental
//      duplications, RepeatMasker)
//   5. Format back to svtk (GatkToSvtkVcf / FORMAT_GATK_VCF_FOR_SVTK,
//      reused from vcfs_cluster_svcluster)
//
// Deliberately not implemented: GATK-SV's cross-*batch* SR evidence flag
// reconciliation (ExtractSRVariantLists/CombineSRBothsidePass/
// SetSRVariantFlags) is skipped -- it exists to combine BOTHSIDES_SUPPORT/
// HIGH_SR_BACKGROUND flags across multiple *batches* of the same caller
// category (GATK-SV's own cohort-mode "batch" concept), which has no
// analog here: we have one VCF per *caller*, not per *batch*, and this
// pipeline doesn't have a first-class batch concept at all (see
// README.md#where-to-simplify-relative-to-upstream-gatk-sv). Also not
// scattered by contig -- same simplification as vcfs_cluster_svcluster,
// see that subworkflow.
//

// Aliased (not plain PLOIDY_TABLE_FROM_PED): this call needs
// retain_female_chr_y=true via conf/modules.config's ext, which is
// selected by process *name* -- aliasing gives it a distinct name so
// vcfs_cluster_svcluster's own PLOIDY_TABLE_FROM_PED call (default false)
// isn't affected by the same config rule.
include { PLOIDY_TABLE_FROM_PED as PLOIDY_TABLE_FROM_PED_COMBINE_BATCHES } from '../../../modules/local/ploidy_table_from_ped/main'
// Both aliased, not plain GATK4_SVCLUSTER: each needs its own overlap
// thresholds via conf/modules.config's ext.args (JOIN: near-zero
// thresholds, exact-duplicate-only; SITES: realistic thresholds, real
// clustering), and vcfs_cluster_svcluster's own GATK4_SVCLUSTER call
// (stage 1, no explicit tuning) must not pick up either rule.
include { GATK4_SVCLUSTER as GATK4_SVCLUSTER_JOIN   } from '../../../modules/nf-core/gatk4/svcluster/main'
include { GATK4_SVCLUSTER as GATK4_SVCLUSTER_SITES  } from '../../../modules/nf-core/gatk4/svcluster/main'
include { GATK4_GROUPEDSVCLUSTER as GATK4_GROUPEDSVCLUSTER_PART1 } from '../../../modules/local/gatk4/grouped_sv_cluster/main'
include { GATK4_GROUPEDSVCLUSTER as GATK4_GROUPEDSVCLUSTER_PART2 } from '../../../modules/local/gatk4/grouped_sv_cluster/main'
include { FORMAT_GATK_VCF_FOR_SVTK } from '../../../modules/local/format_gatk_vcf_for_svtk/main'

workflow VCFS_COMBINE_BATCHES {
    take:
    caller_vcfs                // channel: [ meta, vcf, vcf_tbi ] -- one per caller, meta.id == caller (e.g. from cluster_manta.nf / cluster_wham.nf's clustered_vcf output)
    ped                        // channel: [ meta, ped ]
    contig_list                // channel: [ meta, contig_list ]
    fasta                      // channel: [ meta, fasta ]
    fasta_fai                  // channel: [ meta, fasta_fai ]
    dict                       // channel: [ meta, dict ]
    clustering_config_part1    // channel: [ meta, tsv ]
    stratification_config_part1 // channel: [ meta, tsv ]
    clustering_config_part2    // channel: [ meta, tsv ]
    stratification_config_part2 // channel: [ meta, tsv ]
    track_bed_files             // channel: [ meta, [bed, bed, bed] ] -- one list, e.g. [SR, SD, RM] tracks
    track_names                 // val: List<String>, names matching track_bed_files order (e.g. ['SR','SD','RM'])
    cohort_name                 // val: String, used as the output prefix
    ploidy_script                    // val: absolute path to bin/ploidy_table_from_ped.py
    format_gatk_vcf_for_svtk_script  // val: absolute path to bin/format_gatk_vcf_for_svtk.py

    main:
    // GATK-SV's CombineBatches builds its own ploidy table
    // (wdl/CombineBatches.wdl:110, retain_female_chr_y=true), separate
    // from stage 1's (vcfs_cluster_svcluster's, retain_female_chr_y left
    // at its default false) -- see modules/local/ploidy_table_from_ped
    // for what that flag actually does (rewrites 0 ploidy -> 1, which for
    // females only affects chrY).
    PLOIDY_TABLE_FROM_PED_COMBINE_BATCHES(ped, contig_list.map { meta, c -> c }.collect(), file(ploidy_script))
    ploidy_table = PLOIDY_TABLE_FROM_PED_COMBINE_BATCHES.out.ploidy_table.map { meta, t -> t }.collect()

    // --- Step 1: naive join across callers ---
    // Near-zero overlap thresholds: only merges exact-position
    // duplicates, doesn't attempt real clustering yet. GATK-SV builds
    // this input as flatten([pesr batches, depth batches]); we build it
    // from however many caller_vcfs arrive (currently Manta + Wham).
    // Each stage gets a distinct meta.id (not just cohort_name repeated):
    // GATK4_SVCLUSTER's output filename is derived from meta.id
    // (prefix = task.ext.prefix ?: "${meta.id}"), so if every stage in
    // this chain reused the same cohort_name-only id, each stage's output
    // would collide with the *next* stage's staged input filename (both
    // "cohort.vcf.gz") -- a real bug this way even under -stub-run, found
    // via "Missing output file(s) *.vcf.gz" where the actual problem was
    // the stub overwriting its own staged input in place.
    join_input = caller_vcfs
        .map { meta, vcf, tbi -> [ vcf, tbi ] }
        .collect(flat: false)
        .map { pairs -> [ [id: "${cohort_name}.combine_batches.join_vcfs"], pairs.collect { it[0] }, pairs.collect { it[1] } ] }

    GATK4_SVCLUSTER_JOIN(
        join_input,
        ploidy_table,
        fasta.map { meta, f -> f }.collect(),
        fasta_fai.map { meta, f -> f }.collect(),
        dict.map { meta, d -> d }.collect()
    )

    // --- Step 2: real clustering pass ---
    cluster_input = GATK4_SVCLUSTER_JOIN.out.clustered_vcf
        .join(GATK4_SVCLUSTER_JOIN.out.clustered_vcf_index)
        .map { meta, vcf, tbi -> [ [id: "${cohort_name}.combine_batches.cluster_sites"], [ vcf ], [ tbi ] ] }

    GATK4_SVCLUSTER_SITES(
        cluster_input,
        ploidy_table,
        fasta.map { meta, f -> f }.collect(),
        fasta_fai.map { meta, f -> f }.collect(),
        dict.map { meta, d -> d }.collect()
    )

    // --- Steps 3-4: context-aware re-clustering ---
    // Same per-stage meta.id renaming as steps 1-2, and for the same
    // reason: GATK4_GROUPEDSVCLUSTER's output filename is also derived
    // from meta.id (see modules/local/gatk4/grouped_sv_cluster/main.nf).
    part1_input = GATK4_SVCLUSTER_SITES.out.clustered_vcf
        .join(GATK4_SVCLUSTER_SITES.out.clustered_vcf_index)
        .map { meta, vcf, tbi -> [ [id: "${cohort_name}.combine_batches.recluster_part_1"], vcf, tbi ] }

    GATK4_GROUPEDSVCLUSTER_PART1(
        part1_input,
        ploidy_table,
        fasta.map { meta, f -> f }.collect(),
        fasta_fai.map { meta, f -> f }.collect(),
        dict.map { meta, d -> d }.collect(),
        clustering_config_part1.map { meta, c -> c }.collect(),
        stratification_config_part1.map { meta, c -> c }.collect(),
        track_bed_files.map { meta, beds -> beds }.collect(),
        track_names
    )

    part2_input = GATK4_GROUPEDSVCLUSTER_PART1.out.vcf
        .join(GATK4_GROUPEDSVCLUSTER_PART1.out.vcf_index)
        .map { meta, vcf, tbi -> [ [id: "${cohort_name}.combine_batches.recluster_part_2"], vcf, tbi ] }

    GATK4_GROUPEDSVCLUSTER_PART2(
        part2_input,
        ploidy_table,
        fasta.map { meta, f -> f }.collect(),
        fasta_fai.map { meta, f -> f }.collect(),
        dict.map { meta, d -> d }.collect(),
        clustering_config_part2.map { meta, c -> c }.collect(),
        stratification_config_part2.map { meta, c -> c }.collect(),
        track_bed_files.map { meta, beds -> beds }.collect(),
        track_names
    )

    // --- Step 5: format back to svtk ---
    // GATK-SV uses source="depth" here specifically to match legacy
    // headers (wdl/CombineBatches.wdl's comment: "Use 'depth' as source
    // to match legacy headers"), not because this VCF is depth-only --
    // it's the merged cross-caller cohort site list at this point.
    // Renamed back to a plain cohort_name id for the final output --
    // meta.id up to this point has been a long internal stage-name
    // string, not something we want as the published file's identifier.
    format_input = GATK4_GROUPEDSVCLUSTER_PART2.out.vcf
        .map { meta, vcf -> [ [id: cohort_name], vcf ] }

    FORMAT_GATK_VCF_FOR_SVTK(
        format_input,
        contig_list.map { meta, c -> c }.collect(),
        'depth',
        file(format_gatk_vcf_for_svtk_script)
    )

    combined_vcf = FORMAT_GATK_VCF_FOR_SVTK.out.vcf

    emit:
    combined_vcf  // channel: [ meta, vcf, vcf_tbi ]
}
