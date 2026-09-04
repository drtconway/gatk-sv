//
// Per-variant evidence metrics for FilterBatchSites/AdjudicateSV (svtk
// adjudicate), the random-forest cutoff derivation genotyping needs
// (rf_cutoffs) -- see README.md#path-to-genotyping-the-real-prerequisite-chain
// for why this exists. Mirrors GATK-SV's GenerateBatchMetrics workflow
// (wdl/GenerateBatchMetrics.wdl).
//
// Takes each caller's already-clustered site VCF (vcfs_cluster_svcluster's
// own output, stage 1 -- NOT vcfs_combine_batches' cross-caller merged
// output, confirmed against wdl/GATKSVPipelinePhase1.wdl's own wiring:
// GenerateBatchMetrics.{depth,manta,wham,...}_vcf = ClusterBatch's
// clustered_*_vcf, one call site earlier than CombineBatches). This
// subworkflow does its own simple concat (CONCAT_INPUT_VCFS,
// allow_overlaps) to get one VCF to scatter/annotate -- a different,
// much simpler merge than vcfs_combine_batches' full SVCluster/
// GroupedSVCluster clustering pipeline, which this subworkflow does not
// need or invoke.
//
// Chain: concat all-caller VCFs -> scatter into shards -> per shard
// (FormatVcf svtk->GATK conversion -> SVRegionOverlap segdup/repeat-mask
// annotation -> AggregateSVEvidence PE/SR/BAF-based metrics ->
// AggregateDepthEvidence RD-based metrics) -> concat shards back ->
// AggregateTests (aggregate.py) produces the final evidence_metrics.tsv.
//
// Step 0 (FormatVcf, re-adding ECN) is necessary for the same reason
// vcfs_combine_batches needs it: vcfs_cluster_svcluster's own published
// output has already been converted back to svtk format (ECN stripped),
// so every per-shard GATK walker below needs that conversion re-run
// first -- see this pipeline's CLAUDE.md's own note on the
// svtk-format-vs-GATK-format round trip.
//

include { PLOIDY_TABLE_FROM_PED as PLOIDY_TABLE_FROM_PED_BATCH_METRICS } from '../../../modules/local/ploidy_table_from_ped/main'
include { CONCAT_VCFS as CONCAT_INPUT_VCFS  } from '../../../modules/local/concat_vcfs/main'
include { CONCAT_VCFS as CONCAT_OUTPUT_VCFS } from '../../../modules/local/concat_vcfs/main'
include { SCATTER_VCF } from '../../../modules/local/scatter_vcf/main'
include { TABIX as TABIX_SHARDS } from '../../../modules/local/tabix/main'
include { FORMAT_SVTK_VCF_FOR_GATK as FORMAT_SVTK_VCF_FOR_GATK_BATCH_METRICS } from '../../../modules/local/format_svtk_vcf_for_gatk/main'
include { TABIX as TABIX_FORMATTED_SHARDS } from '../../../modules/local/tabix/main'
include { GATK_SVREGIONOVERLAP } from '../../../modules/local/gatk4/sv_region_overlap/main'
include { GATK_AGGREGATESVEVIDENCE } from '../../../modules/local/gatk4/aggregate_sv_evidence/main'
include { GATK_AGGREGATEDEPTHEVIDENCE } from '../../../modules/local/gatk4/aggregate_depth_evidence/main'
include { AGGREGATE_TESTS } from '../../../modules/local/aggregate_tests/main'

workflow VCFS_GENERATE_BATCH_METRICS {
    take:
    caller_vcfs        // channel: [ meta, vcf, vcf_tbi ] -- one per caller, meta.id == caller (e.g. from vcfs_cluster_svcluster's clustered_vcf output)
    ped                // channel: [ meta, ped ]
    contig_list        // channel: [ meta, contig_list ]
    dict               // channel: [ meta, dict ]
    pe_file            // channel: [ meta, pe_txt_gz, pe_txt_gz_tbi ] -- panel-wide, from vcfs_merge_evidence
    sr_file            // channel: [ meta, sr_txt_gz, sr_txt_gz_tbi ] -- panel-wide
    baf_file           // channel: [ meta, baf_txt_gz, baf_txt_gz_tbi ] -- panel-wide
    rd_file            // channel: [ meta, rd_txt_gz, rd_txt_gz_tbi ] -- panel-wide, from vcfs_merge_read_counts
    median_coverage    // channel: [ meta, median_cov_bed ]
    segdups            // channel: [ meta, bed, bed_tbi ]
    rmsk               // channel: [ meta, bed, bed_tbi ]
    chr_x              // val: String, e.g. 'chrX'
    chr_y              // val: String, e.g. 'chrY'
    cohort_name        // val: String, used as the output prefix
    records_per_shard   // val: Integer
    ploidy_script                    // val: absolute path to bin/ploidy_table_from_ped.py
    format_svtk_vcf_for_gatk_script  // val: absolute path to bin/format_svtk_vcf_for_gatk.py
    aggregate_script                 // val: absolute path to bin/aggregate.py

    main:
    // GenerateBatchMetrics builds its own ploidy table
    // (wdl/GenerateBatchMetrics.wdl:74, retain_female_chr_y left at
    // its default false) -- separate from vcfs_combine_batches' own
    // (retain_female_chr_y=true), since these are two different WDL call
    // sites, not the same table reused.
    PLOIDY_TABLE_FROM_PED_BATCH_METRICS(ped, contig_list.map { meta, c -> c }.collect(), file(ploidy_script))
    ploidy_table = PLOIDY_TABLE_FROM_PED_BATCH_METRICS.out.ploidy_table.map { meta, t -> t }.collect()

    // --- Concat all callers' clustered VCFs into one, then scatter into shards ---
    concat_input = caller_vcfs
        .map { meta, vcf, tbi -> [ vcf, tbi ] }
        .collect(flat: false)
        .map { pairs -> [ [id: "${cohort_name}.batch_metrics.concat_input_vcfs"], pairs.collect { it[0] }, pairs.collect { it[1] } ] }

    CONCAT_INPUT_VCFS(concat_input)

    SCATTER_VCF(CONCAT_INPUT_VCFS.out.vcf, records_per_shard)

    // Each shard needs its own tuple [meta, vcf] (no index -- FormatVcf's
    // own module takes an un-indexed VCF, same as vcfs_combine_batches'
    // equivalent call), with a distinct per-shard id (not just
    // cohort_name repeated) so this scatter's per-process outputs don't
    // collide with each other -- same "process-name-derived filenames
    // collide across a multi-stage chain" pitfall as
    // vcfs_combine_batches (see that subworkflow's own note). shard's own
    // filename (from SCATTER_VCF, meta.id "${cohort_name}.batch_metrics.
    // concat_input_vcfs") already carries a unique-enough name per shard
    // (via its own "shard_NNNNNN" suffix) -- reuse the shard's own
    // baseName directly as the new meta.id rather than prepending
    // cohort_name/"batch_metrics" a second time.
    shards = SCATTER_VCF.out.shards
        .flatMap { meta, shard_list ->
            def list = shard_list instanceof List ? shard_list : [ shard_list ]
            list.collect { shard -> [ [id: shard.baseName - '.vcf'], shard ] }
        }

    // SCATTER_VCF's own bcftools +scatter doesn't index its shards --
    // index explicitly, since FORMAT_SVTK_VCF_FOR_GATK's underlying
    // script needs a queryable (indexed) VCF -- found for real (not
    // caught by any earlier per-module test in isolation, since this
    // subworkflow is the first place a SCATTER_VCF shard feeds directly
    // into FORMAT_SVTK_VCF_FOR_GATK).
    TABIX_SHARDS(shards)

    // --- Per shard: FormatVcf -> SVRegionOverlap -> AggregateSVEvidence -> AggregateDepthEvidence ---
    FORMAT_SVTK_VCF_FOR_GATK_BATCH_METRICS(
        TABIX_SHARDS.out.vcf.map { meta, vcf, tbi -> [ meta, vcf ] },
        ploidy_table.map { t -> [ [id: 'ploidy'], t ] },
        file(format_svtk_vcf_for_gatk_script)
    )

    // FORMAT_SVTK_VCF_FOR_GATK doesn't emit a .tbi -- index explicitly,
    // same as vcfs_combine_batches' own equivalent step.
    TABIX_FORMATTED_SHARDS(FORMAT_SVTK_VCF_FOR_GATK_BATCH_METRICS.out.vcf)
    formatted_vcf = TABIX_FORMATTED_SHARDS.out.vcf

    // segdups/rmsk are single, static tracks broadcast across every
    // shard below (formatted_vcf is scattered, one call per shard).
    // .collect() flattens a list-valued channel item just like .combine()
    // does (confirmed empirically, not assumed -- see
    // vcfs_merge_read_counts's own note on the analogous .combine() case):
    // collecting a channel of single [bed, tbi] pairs gives [bed, tbi]
    // (size 2), not [[bed, tbi]] (size 1). Wrapping in an extra list
    // before collecting keeps it as one distinguishable tuple element
    // afterward, then unwrap the one-element outer list back to the bare
    // [bed, tbi] pair GATK_SVREGIONOVERLAP's tuple input expects.
    segdups_value = segdups.map { meta, bed, tbi -> [ [ bed, tbi ] ] }.collect().map { it[0] }
    rmsk_value = rmsk.map { meta, bed, tbi -> [ [ bed, tbi ] ] }.collect().map { it[0] }

    GATK_SVREGIONOVERLAP(
        formatted_vcf,
        dict.map { meta, d -> d }.collect(),
        segdups_value,
        rmsk_value
    )

    // Each stage below gets a distinct meta.id suffix, not the bare shard
    // id reused unchanged: GATK_SVREGIONOVERLAP/GATK_AGGREGATESVEVIDENCE/
    // GATK_AGGREGATEDEPTHEVIDENCE all derive their own output filename
    // from meta.id (task.ext.prefix ?: "${meta.id}"), so reusing the same
    // id across this whole per-shard chain makes each stage's output
    // filename identical to the *next* stage's staged input filename --
    // found for real (not caught by any single-module test in isolation):
    // AggregateSVEvidence's own -V/-O both resolved to the same path,
    // silently processing 0 variants rather than erroring, since its
    // input and output happened to be the same file. Same root cause as
    // vcfs_combine_batches' own documented instance of this pitfall.
    svregionoverlap_vcf = GATK_SVREGIONOVERLAP.out.vcf
        .map { meta, vcf, tbi -> [ [id: "${meta.id}.svregionoverlap"], vcf, tbi ] }

    // pe_file/sr_file/baf_file/rd_file are each a single panel-wide value
    // (one emission total, from vcfs_merge_evidence/vcfs_merge_read_counts
    // upstream) that must broadcast across every shard the same way
    // dict/segdups/rmsk do above -- GATK_SVREGIONOVERLAP.out.vcf is
    // scattered, one emission per shard. Same wrap-collect-unwrap as
    // segdups/rmsk above, preserving the full [meta, file, index] tuple
    // shape each of these modules' own tuple inputs expects (unlike
    // segdups/rmsk, whose modules take a bare [file, index] pair with no
    // meta).
    pe_value  = pe_file.map { meta, f, tbi -> [ [ meta, f, tbi ] ] }.collect().map { it[0] }
    sr_value  = sr_file.map { meta, f, tbi -> [ [ meta, f, tbi ] ] }.collect().map { it[0] }
    baf_value = baf_file.map { meta, f, tbi -> [ [ meta, f, tbi ] ] }.collect().map { it[0] }
    rd_value  = rd_file.map { meta, f, tbi -> [ [ meta, f, tbi ] ] }.collect().map { it[0] }

    GATK_AGGREGATESVEVIDENCE(
        svregionoverlap_vcf,
        median_coverage.map { meta, m -> m }.collect(),
        ploidy_table.map { t -> [ [id: 'ploidy'], t ] },
        pe_value,
        sr_value,
        baf_value,
        chr_x,
        chr_y
    )

    // Same per-stage rename as svregionoverlap_vcf above, same reason.
    aggregatesvevidence_vcf = GATK_AGGREGATESVEVIDENCE.out.vcf
        .map { meta, vcf, tbi -> [ [id: "${meta.id}.aggregatesvevidence"], vcf, tbi ] }

    GATK_AGGREGATEDEPTHEVIDENCE(
        aggregatesvevidence_vcf,
        median_coverage.map { meta, m -> m }.collect(),
        rd_value
    )

    // --- Concat shards back, then flatten to the final metrics.tsv ---
    // Sort by filename (not full path!) before concatenating: unlike
    // WDL/Cromwell's scatter (whose output array is always in
    // scatter-index order), Nextflow's .collect() gathers a scattered
    // channel's emissions in parallel-completion order, not shard order
    // -- found for real (not assumed): CONCAT_OUTPUT_VCFS failed on real
    // shard output with "Unsorted positions ... 400 followed by 100"
    // because shard_000001 happened to finish before shard_000000. A
    // first attempt at this fix sorted on it[0].toString() (the full
    // absolute path, e.g. ".../work/c0/c41846.../shard_000000.vcf.gz")
    // and silently did nothing useful -- it sorts by the *work directory
    // hash prefix* (random per run), not the shard filename suffix,
    // confirmed empirically after the first fix attempt didn't change
    // the output order at all. .name (just the filename) is what
    // actually needs sorting, since SCATTER_VCF's own shard names are
    // zero-padded ("shard_000000", "shard_000001", ...) and lexically
    // sortable. Same fix wdl/TasksMakeCohortVcf.wdl's own ConcatVcfs task
    // offers via its sort_vcf_list option (not set at this WDL's own
    // ConcatOutputVcfs call site, since Cromwell's scatter output array
    // doesn't need it).
    concat_output = GATK_AGGREGATEDEPTHEVIDENCE.out.vcf
        .map { meta, vcf, tbi -> [ vcf, tbi ] }
        .collect(flat: false)
        .map { pairs -> pairs.sort { it[0].name } }
        .map { pairs -> [ [id: "${cohort_name}.batch_metrics.concat_output_vcfs"], pairs.collect { it[0] }, pairs.collect { it[1] } ] }

    CONCAT_OUTPUT_VCFS(concat_output)

    AGGREGATE_TESTS(
        CONCAT_OUTPUT_VCFS.out.vcf,
        file(aggregate_script)
    )

    metrics = AGGREGATE_TESTS.out.metrics

    emit:
    metrics       // channel: [ meta, metrics_tsv ]
    ploidy_table  // channel: [ ploidy_tsv ] -- single-value list, GenerateBatchMetrics.wdl also emits this as an output for downstream reuse
}
