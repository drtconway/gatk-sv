//
// Annotate a VCF's variants with RD (read-depth)-based support metrics,
// for AggregateTests downstream. Wraps GATK's AggregateDepthEvidence
// walker -- no nf-core module exists for it. Reproduces
// wdl/GenerateBatchMetrics.wdl's AggregateDepthEvidence task.
//
// Same GATK-SV-fork-only container as modules/local/gatk4/
// sv_region_overlap -- see that module's own note on why this differs
// from every other gatk4/* module here (gatk4-main_gcnvkernel doesn't
// have this walker at all).
//
process GATK_AGGREGATEDEPTHEVIDENCE {
    tag "${meta.id}"
    label 'process_single'

    container 'us.gcr.io/broad-dsde-methods/gatk-sv/gatk:mw-gatk-sv-53d5c2d'

    input:
    tuple val(meta), path(vcf), path(vcf_index)
    path(median_coverage)
    tuple val(rd_meta), path(rd_file), path(rd_file_index)

    output:
    tuple val(meta), path("*.vcf.gz"), path("*.vcf.gz.tbi"), emit: vcf

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    def avail_mem = 3072
    if (!task.memory) {
        log.info('[GATK AggregateDepthEvidence] Available memory not known - defaulting to 3GB. Specify process memory requirements to change this.')
    }
    else {
        avail_mem = (task.memory.mega * 0.8).intValue()
    }
    """
    gatk --java-options "-Xmx${avail_mem}M -XX:-UsePerfData" \\
        AggregateDepthEvidence \\
        -V ${vcf} \\
        -O ${prefix}.vcf.gz \\
        --median-coverage ${median_coverage} \\
        --rd-file ${rd_file} \\
        ${args}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo "" | gzip > ${prefix}.vcf.gz
    touch ${prefix}.vcf.gz.tbi
    """
}
