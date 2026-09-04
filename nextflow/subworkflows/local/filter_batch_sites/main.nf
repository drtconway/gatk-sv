//
// RF-score-based sites filtering: takes each caller's own sites VCF (from
// vcfs_cluster_svcluster, stage 1) plus AdjudicateSV's scores/cutoffs and
// GenerateBatchMetrics' own evidence_metrics.tsv, filters each caller's
// VCF to score-passing records, then merges across callers into one VCF
// -- GenotypeBatch's real `vcf` input upstream (confirmed against
// wdl/GATKSVPipelineBatch.wdl's own wiring: GenotypeBatch.vcf =
// MergePesrDepthVcfs.concat_vcf, itself select_all([filtered_pesr_vcf,
// filtered_depth_vcf]) -- FilterBatchSites' own FilterAnnotateVcf output,
// not CombineBatches' or GenerateBatchMetrics' raw sites VCF). Mirrors
// GATK-SV's FilterBatchSites workflow (wdl/FilterBatchSites.wdl), minus
// its own PlotSVCountsPerSample (outlier-sample QC plotting -- not on
// the genotyping path, same scope decision as AdjudicateSV's own module,
// see README.md's own note).
//
// This pipeline only has Manta/Wham callers so far (both PESR-type; no
// depth caller yet -- cn.MOPS/gCNV are both "Not yet started", see
// README.md), so the merge step here is simpler than upstream's
// Array[File?] depth_vcf/dragen_vcf/manta_vcf/melt_vcf/scramble_vcf/
// wham_vcf scatter: takes however many caller_vcfs arrive (currently
// Manta + Wham), same "build from whatever's actually wired up" pattern
// vcfs_combine_batches/vcfs_generate_batch_metrics already use for their
// own caller_vcfs inputs.
//

include { FILTER_ANNOTATE_VCF } from '../../../modules/local/filter_annotate_vcf/main'
include { TABIX as TABIX_FILTERED_VCFS } from '../../../modules/local/tabix/main'
include { CONCAT_VCFS as MERGE_PESR_DEPTH_VCFS } from '../../../modules/local/concat_vcfs/main'

workflow FILTER_BATCH_SITES {
    take:
    caller_vcfs   // channel: [ meta, vcf, vcf_tbi ] -- one per caller, meta.id == caller (e.g. from vcfs_cluster_svcluster's clustered_vcf output)
    metrics       // channel: [ meta, metrics_tsv ] -- from vcfs_generate_batch_metrics/AggregateTests
    scores        // channel: [ meta, scores_tsv ] -- from AdjudicateSV
    cutoffs       // channel: [ meta, cutoffs_tsv ] -- from AdjudicateSV (this is rf_cutoffs)
    cohort_name   // val: String, used as the output prefix
    rewrite_sr_coords_script     // val: absolute path to bin/rewrite_SR_coords.py
    annotate_rf_evidence_script  // val: absolute path to bin/annotate_RF_evidence.py

    main:
    // metrics/scores/cutoffs are each a single panel-wide value broadcast
    // across every caller below -- same wrap-collect-unwrap pattern as
    // vcfs_generate_batch_metrics/genotype_batch (see those subworkflows'
    // own notes on why .collect() alone flattens a multi-element tuple).
    metrics_value = metrics.map { meta, m -> [ [ meta, m ] ] }.collect().map { it[0] }
    scores_value  = scores.map { meta, s -> [ [ meta, s ] ] }.collect().map { it[0] }
    cutoffs_value = cutoffs.map { meta, c -> [ [ meta, c ] ] }.collect().map { it[0] }

    FILTER_ANNOTATE_VCF(
        caller_vcfs,
        metrics_value,
        scores_value,
        cutoffs_value,
        file(rewrite_sr_coords_script),
        file(annotate_rf_evidence_script)
    )

    // FILTER_ANNOTATE_VCF doesn't emit a .tbi -- index explicitly, same
    // pattern as vcfs_combine_batches'/vcfs_generate_batch_metrics' own
    // equivalent steps (FORMAT_SVTK_VCF_FOR_GATK doesn't emit one
    // either). wdl/GATKSVPipelineBatch.wdl's own MergePesrDepthVcfs call
    // builds an index array the same way (merge_vcfs_[i] + ".tbi"),
    // since FilterAnnotateVcf's real output has none either -- but a
    // bare file("${vcf}.tbi") reference here would just be a guess that
    // the file exists, not a real index; found and fixed this exact
    // pattern earlier in vcfs_combine_batches/vcfs_generate_batch_metrics
    // and it's the same mistake to avoid repeating.
    TABIX_FILTERED_VCFS(FILTER_ANNOTATE_VCF.out.annotated_vcf)

    merge_input = TABIX_FILTERED_VCFS.out.vcf
        .map { meta, vcf, tbi -> [ vcf, tbi ] }
        .collect(flat: false)
        .map { pairs -> [ [id: "${cohort_name}.merge_pesr_depth"], pairs.collect { it[0] }, pairs.collect { it[1] } ] }

    MERGE_PESR_DEPTH_VCFS(merge_input)

    emit:
    merged_vcf = MERGE_PESR_DEPTH_VCFS.out.vcf  // channel: [ meta, vcf, vcf_tbi ] -- GenotypeBatch's own `vcf` input
}
