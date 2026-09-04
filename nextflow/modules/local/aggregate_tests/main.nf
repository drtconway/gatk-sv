//
// Flatten a VCF's per-variant evidence annotations (from
// AggregateSVEvidence/AggregateDepthEvidence/SVRegionOverlap upstream)
// into the flat evidence_metrics.tsv that FilterBatchSites/AdjudicateSV
// (svtk adjudicate) needs. Wraps GATK-SV's own aggregate.py (vendored in
// ../../../bin/, adapted -- see bin/README.md's own note on why it no
// longer imports svtk.utils -- see wdl/GenerateBatchMetrics.wdl's
// AggregateTests task, which this reproduces).
//
// Deliberately not implemented: the WDL task's optional
// outlier_sample_ids input (excludes given samples' variants from
// training) -- no caller in this pipeline has an outlier-sample list to
// pass yet (FilterBatchSamples' outlier detection isn't built here).
// Revisit if/when that exists rather than threading an unused optional
// input through now.
//
// The script is passed in as an explicit `path` input rather than relying
// on Nextflow's bin/-auto-PATH mechanism -- see
// modules/local/ploidy_table_from_ped/main.nf's own note on why that
// mechanism never fires for any workflows/*.nf entry point here.
//
process AGGREGATE_TESTS {
    tag "${meta.id}"
    label 'process_single'

    container 'drtomc/gatk-sv-nf-aggregate-tests:0.1.0'

    input:
    tuple val(meta), path(vcf), path(vcf_index)
    path(script)

    output:
    tuple val(meta), path("*.metrics.tsv"), emit: metrics

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    python3 ${script} \\
        -v ${vcf} \\
        ${prefix}.metrics.tsv
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.metrics.tsv
    """
}
