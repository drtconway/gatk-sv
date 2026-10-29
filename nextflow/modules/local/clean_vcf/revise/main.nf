//
// CleanVcf's three GATK-SV-fork walkers for one contig, chained:
// SVReviseOverlappingCnvs, SVReviseMultiallelicCnvs,
// SVReviseOverlappingMultiallelics (copy-number consistency between
// overlapping CNVs, multiallelic CNV calls). Reproduces
// wdl/CleanVcfChromosome.wdl's three CleanVcfRevise* tasks in one. Same
// GATK-SV-fork image as GATK_TRAINSVGENOTYPING; the sv-pipeline image's
// GATK jar doesn't have these walkers.
//
// No outlier-sample list (upstream's optional --outlier-samples): this
// pipeline has no FilterBatchSamples yet.
//
process CLEANVCF_REVISE {
    tag "${meta.id}"
    label 'process_single'

    container 'us.gcr.io/broad-dsde-methods/gatk-sv/gatk:mw-gatk-sv-53d5c2d'

    input:
    tuple val(meta), path(vcf), path(vcf_index)

    output:
    tuple val(meta), path("${task.ext.prefix ?: meta.id}.vcf.gz"), path("${task.ext.prefix ?: meta.id}.vcf.gz.tbi"), emit: vcf

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def avail_mem = task.memory ? (task.memory.mega * 0.8).intValue() : 3072
    """
    gatk --java-options "-Xmx${avail_mem}m" SVReviseOverlappingCnvs \\
        -V ${vcf} \\
        -O revise_overlapping_cnvs.vcf.gz

    gatk --java-options "-Xmx${avail_mem}m" SVReviseMultiallelicCnvs \\
        -V revise_overlapping_cnvs.vcf.gz \\
        -O revise_multiallelic_cnvs.vcf.gz

    gatk --java-options "-Xmx${avail_mem}m" SVReviseOverlappingMultiallelics \\
        -V revise_multiallelic_cnvs.vcf.gz \\
        -O ${prefix}.vcf.gz
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo "" | gzip > ${prefix}.vcf.gz
    touch ${prefix}.vcf.gz.tbi
    """
}
