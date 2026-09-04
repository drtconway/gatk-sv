//
// Random-forest adjudication of per-variant evidence metrics
// (GenerateBatchMetrics/AggregateTests' own evidence_metrics.tsv) into
// scores and rf_cutoffs -- the direct input TrainSVGenotyping/GenotypeSVs
// need for genotyping. Wraps GATK-SV's own svtk.adjudicate package
// (vendored in ../../../bin/, adapted -- see bin/README.md's own note on
// why this calls the vendored logic directly rather than installing
// svtk itself). Reproduces wdl/FilterBatchSites.wdl's AdjudicateSV task.
//
// Deliberately not implemented: FilterBatchSites' own FilterAnnotateVcf
// (per-caller VCF filtering using these cutoffs/scores) and
// PlotSVCountsPerSample (outlier-sample QC plotting) -- neither is on
// the path to genotyping, which only needs this task's own cutoffs
// output (rf_cutoffs). Revisit together if/when this pipeline needs
// per-caller sites-filtered VCFs or outlier-sample QC for its own sake.
//
// The scripts are passed in as explicit `path` inputs rather than
// relying on Nextflow's bin/-auto-PATH mechanism -- see
// modules/local/ploidy_table_from_ped/main.nf's own note on why that
// mechanism never fires for any workflows/*.nf entry point here.
//
process ADJUDICATE_SV {
    tag "${meta.id}"
    label 'process_single'

    container 'drtomc/gatk-sv-nf-adjudicate:0.1.0'

    input:
    tuple val(meta), path(metrics)
    path(adjudicate_script)
    path(adjudicate_sv_script)
    path(random_forest_script)
    path(labelers_script)

    output:
    tuple val(meta), path("*.scores"), emit: scores
    tuple val(meta), path("*.cutoffs"), emit: cutoffs
    tuple val(meta), path("*.RF_intermediate_files.tar.gz"), emit: rf_intermediate_files

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    python3 ${adjudicate_script} ${metrics} ${prefix}.scores ${prefix}.cutoffs

    mkdir ${prefix}.RF_intermediate_files
    mv *_trainable.txt ${prefix}.RF_intermediate_files/
    mv *_testable.txt ${prefix}.RF_intermediate_files/
    tar -czf ${prefix}.RF_intermediate_files.tar.gz ${prefix}.RF_intermediate_files
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.scores ${prefix}.cutoffs
    mkdir ${prefix}.RF_intermediate_files
    tar -czf ${prefix}.RF_intermediate_files.tar.gz ${prefix}.RF_intermediate_files
    """
}
