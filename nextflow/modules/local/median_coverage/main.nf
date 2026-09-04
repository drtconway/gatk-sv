//
// Compute each sample's median bin coverage from the panel-wide bincov
// matrix (vcfs_merge_read_counts's output), for GATK's
// AggregateSVEvidence/AggregateDepthEvidence and TrainSVGenotyping/
// GenotypeSVs (all take a --median-coverage file with this exact shape).
// Wraps GATK-SV's own medianCoverage.R (vendored in ../../../bin/, see
// wdl/MedianCov.wdl's CalcMedCov task, which this reproduces).
//
// The script is passed in as an explicit `path` input rather than relying
// on Nextflow's bin/-auto-PATH mechanism -- see
// modules/local/ploidy_table_from_ped/main.nf's own note on why that
// mechanism never fires for any workflows/*.nf entry point here.
//
process MEDIAN_COVERAGE {
    tag "${meta.id}"
    label 'process_single'

    container 'drtomc/gatk-sv-nf-median-coverage:0.1.0'

    input:
    tuple val(meta), path(bincov_matrix)
    path(script)

    output:
    tuple val(meta), path("*.medianCov.transposed.bed"), emit: median_coverage

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    # Named distinctly from the input (not e.g. "fixed.bed"): a bincov
    # matrix input that happened to already be named that would have
    # this decompress-redirect truncate its own input before reading it.
    zcat -f ${bincov_matrix} > tmp_fixed.bed

    Rscript ${script} -H tmp_fixed.bed medianCov.bed

    # Reproduces wdl/MedianCov.wdl's own post-processing exactly: read
    # positionally (check.names=FALSE, no header=TRUE) since
    # medianCoverage.R's own output header starts with "#sample_id" --
    # R's default comment.char="#" would otherwise swallow that whole
    # header line as a comment. Columns 1-2 are #sample_id/Med_withZeros;
    # column 3 (Med_withoutZeros) is dropped here too, matching upstream,
    # even though it looks like the more useful of the two -- not
    # something to silently \"fix\" without checking why GATK-SV's own
    # task picked the with-zeros column.
    Rscript -e '
    x <- read.table("medianCov.bed", check.names=FALSE)
    xtransposed <- t(x[, c(1, 2)])
    write.table(xtransposed, file="${prefix}.medianCov.transposed.bed", sep="\\t", row.names=FALSE, col.names=FALSE, quote=FALSE)
    '
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.medianCov.transposed.bed
    """
}
