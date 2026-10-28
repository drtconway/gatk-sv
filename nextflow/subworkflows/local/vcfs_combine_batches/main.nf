//
// Merge genotype_batch's genotyped PESR and depth VCFs into one cohort VCF,
// clustering the same variant called by different callers (Manta and Wham
// records are still separate after genotyping) into one record. Mirrors
// GATK-SV's CombineBatches (wdl/CombineBatches.wdl), the first step of
// MakeCohortVcf, which upstream also runs on GenotypeBatch's output.
//
// Steps, matching CombineBatches:
//   1. Naive join (SVCluster, near-zero overlap thresholds: only merges
//      exact-position duplicates)
//   2. Real clustering pass (SVCluster, realistic overlap thresholds)
//   3+4. Two rounds of context-aware re-clustering (GroupedSVCluster,
//      stratified by simple repeat, segdup and RepeatMasker tracks)
//   5. ExtractSRVariantLists: HIGH_SR_BACKGROUND / BOTHSIDES_SUPPORT ID
//      lists, which ResolveComplexVariants and CleanVcf take as inputs
//   6. Format to svtk (GatkToSvtkVcf: source "depth", PASS set, AC/AF/AN,
//      SR flags and SR1POS/SR2POS dropped, as upstream)
//
// The inputs are already GATK format (GenotypeSVs writes ECN), so unlike
// the earlier version of this subworkflow, no svtk->GATK conversion first.
//
// Not implemented, deliberately: upstream's cross-batch SR flag
// reconciliation (CombineSRBothsidePass, CombineBackgroundFail,
// SetSRVariantFlags). With a single batch it is close to an identity: the
// background-fail list is unchanged, and BOTHSIDES_SUPPORT is only dropped
// from variants with no non-ref genotype, which CleanVcf removes anyway.
// Also not scattered by contig, unlike upstream (same simplification as
// vcfs_cluster_svcluster).
//

// Aliased: retain_female_chr_y=true is set by process name in
// conf/modules.config, and stage 1's PLOIDY_TABLE_FROM_PED must not get it.
include { PLOIDY_TABLE_FROM_PED as PLOIDY_TABLE_FROM_PED_COMBINE_BATCHES } from '../../../modules/local/ploidy_table_from_ped/main'
include { PED_SEXES } from '../../../modules/local/ped_sexes/main'
// Aliased: each gets its own overlap thresholds via ext.args, which stage
// 1's GATK4_SVCLUSTER must not pick up.
include { GATK4_SVCLUSTER as GATK4_SVCLUSTER_JOIN   } from '../../../modules/nf-core/gatk4/svcluster/main'
include { GATK4_SVCLUSTER as GATK4_SVCLUSTER_SITES  } from '../../../modules/nf-core/gatk4/svcluster/main'
include { GATK4_GROUPEDSVCLUSTER as GATK4_GROUPEDSVCLUSTER_PART1 } from '../../../modules/local/gatk4/grouped_sv_cluster/main'
include { GATK4_GROUPEDSVCLUSTER as GATK4_GROUPEDSVCLUSTER_PART2 } from '../../../modules/local/gatk4/grouped_sv_cluster/main'
include { EXTRACT_SR_VARIANT_LISTS } from '../../../modules/local/extract_sr_variant_lists/main'
// Aliased: needs its own ext.args (field removal, --set-pass) and
// publishDir, separate from stage 1's FORMAT_GATK_VCF_FOR_SVTK.
include { FORMAT_GATK_VCF_FOR_SVTK as FORMAT_GATK_VCF_FOR_SVTK_COMBINE_BATCHES } from '../../../modules/local/format_gatk_vcf_for_svtk/main'

workflow VCFS_COMBINE_BATCHES {
    take:
    pesr_vcf                    // channel: [ meta, vcf, vcf_tbi ] -- genotype_batch's genotyped_pesr_vcf (GATK format)
    depth_vcf                   // channel: [ meta, vcf, vcf_tbi ] -- genotype_batch's genotyped_depth_vcf (empty until a depth caller exists)
    ped                         // channel: [ meta, ped ]
    contig_list                 // channel: [ meta, contig_list ]
    fasta                       // channel: [ meta, fasta ]
    fasta_fai                   // channel: [ meta, fasta_fai ]
    dict                        // channel: [ meta, dict ]
    clustering_config_part1     // channel: [ meta, tsv ]
    stratification_config_part1 // channel: [ meta, tsv ]
    clustering_config_part2     // channel: [ meta, tsv ]
    stratification_config_part2 // channel: [ meta, tsv ]
    track_bed_files             // channel: [ meta, [bed, bed, bed] ] -- one list, e.g. [SR, SD, RM] tracks
    track_names                 // val: List<String>, names matching track_bed_files order (e.g. ['SR','SD','RM'])
    cohort_name                 // val: String, used as the output prefix and variant ID prefix
    ploidy_script                    // val: absolute path to bin/ploidy_table_from_ped.py
    format_gatk_vcf_for_svtk_script  // val: absolute path to bin/format_gatk_vcf_for_svtk.py

    main:
    // Sexes only, so pedigree-only PED edits stay cached (see PED_SEXES).
    PED_SEXES(ped)
    PLOIDY_TABLE_FROM_PED_COMBINE_BATCHES(PED_SEXES.out.ped, contig_list.map { meta, c -> c }.collect(), file(ploidy_script))
    ploidy_table = PLOIDY_TABLE_FROM_PED_COMBINE_BATCHES.out.ploidy_table.map { meta, t -> t }.collect()

    fasta_v     = fasta.map { meta, f -> f }.collect()
    fasta_fai_v = fasta_fai.map { meta, f -> f }.collect()
    dict_v      = dict.map { meta, d -> d }.collect()

    // Each stage gets its own meta.id: output filenames derive from it, and
    // a reused id would write through the next stage's staged input.
    // meta.cohort carries the variant ID prefix (GATK4_SVCLUSTER_SITES'
    // ext.args).

    // --- Step 1: naive join, PESR + depth (upstream: flatten([pesr, depth])) ---
    join_input = pesr_vcf.mix(depth_vcf)
        .map { meta, vcf, tbi -> [ vcf, tbi ] }
        .collect(flat: false)
        .map { pairs -> pairs.sort { it[0].name } }
        .map { pairs -> [ [id: "${cohort_name}.combine_batches.join_vcfs", cohort: cohort_name], pairs.collect { it[0] }, pairs.collect { it[1] } ] }

    GATK4_SVCLUSTER_JOIN(join_input, ploidy_table, fasta_v, fasta_fai_v, dict_v)

    // --- Step 2: real clustering pass ---
    cluster_input = GATK4_SVCLUSTER_JOIN.out.clustered_vcf
        .join(GATK4_SVCLUSTER_JOIN.out.clustered_vcf_index)
        .map { meta, vcf, tbi -> [ [id: "${cohort_name}.combine_batches.cluster_sites", cohort: cohort_name], [ vcf ], [ tbi ] ] }

    GATK4_SVCLUSTER_SITES(cluster_input, ploidy_table, fasta_v, fasta_fai_v, dict_v)

    // --- Steps 3-4: context-aware re-clustering ---
    part1_input = GATK4_SVCLUSTER_SITES.out.clustered_vcf
        .join(GATK4_SVCLUSTER_SITES.out.clustered_vcf_index)
        .map { meta, vcf, tbi -> [ [id: "${cohort_name}.combine_batches.recluster_part_1"], vcf, tbi ] }

    GATK4_GROUPEDSVCLUSTER_PART1(
        part1_input,
        ploidy_table, fasta_v, fasta_fai_v, dict_v,
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
        ploidy_table, fasta_v, fasta_fai_v, dict_v,
        clustering_config_part2.map { meta, c -> c }.collect(),
        stratification_config_part2.map { meta, c -> c }.collect(),
        track_bed_files.map { meta, beds -> beds }.collect(),
        track_names
    )

    clustered = GATK4_GROUPEDSVCLUSTER_PART2.out.vcf.join(GATK4_GROUPEDSVCLUSTER_PART2.out.vcf_index)

    // --- Step 5: SR flag lists, from the GATK-format VCF (step 6 drops them) ---
    EXTRACT_SR_VARIANT_LISTS(
        clustered.map { meta, vcf, tbi -> [ [id: "${cohort_name}.combine_batches"], vcf, tbi ] }
    )

    // --- Step 6: format to svtk ---
    // source "depth" only matches upstream's legacy headers; this is the
    // merged PESR + depth VCF.
    FORMAT_GATK_VCF_FOR_SVTK_COMBINE_BATCHES(
        clustered.map { meta, vcf, tbi -> [ [id: "${cohort_name}.combine_batches"], vcf ] },
        contig_list.map { meta, c -> c }.collect(),
        'depth',
        file(format_gatk_vcf_for_svtk_script)
    )

    emit:
    combined_vcf    = FORMAT_GATK_VCF_FOR_SVTK_COMBINE_BATCHES.out.vcf  // channel: [ meta, vcf, vcf_tbi ]
    background_fail = EXTRACT_SR_VARIANT_LISTS.out.background_fail      // channel: [ meta, txt ]
    bothside_pass   = EXTRACT_SR_VARIANT_LISTS.out.bothside_pass        // channel: [ meta, txt ]
}
