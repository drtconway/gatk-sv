//
// One contig's records, and the CNV intervals of its complex SVs (and
// unresolved inversion single-enders) to depth-genotype. Reproduces
// wdl/GenotypeCpxCnvs.wdl's GetCpxCnvIntervals, plus the per-contig split
// upstream gets from ReshardVcf.
//
process CPX_CNV_INTERVALS {
    tag "${meta.id}"
    label 'process_single'

    container 'us.gcr.io/broad-dsde-methods/gatk-sv/sv-pipeline:2026-06-02-v1.1-24f0fd08'

    input:
    tuple val(meta), path(vcf), path(vcf_index)  // meta.contig

    output:
    tuple val(meta), path("${task.ext.prefix ?: meta.id}.vcf.gz"), path("${task.ext.prefix ?: meta.id}.vcf.gz.tbi"), emit: vcf
    tuple val(meta), path("${task.ext.prefix ?: meta.id}.complex_CNV_intervals_to_test.bed.gz"), emit: intervals

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    // The script's own temp dir is mktemp -d (/tmp); its files are small.
    """
    bcftools view --no-version -r ${meta.contig} ${vcf} -Oz -o ${prefix}.vcf.gz
    tabix ${prefix}.vcf.gz

    bash /opt/sv-pipeline/04_variant_resolution/scripts/gather_cpx_intervals_for_rd_gt.sh \\
        ${prefix}.vcf.gz \\
        ${prefix}.complex_CNV_intervals_to_test.bed.gz
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo "" | gzip > ${prefix}.vcf.gz
    touch ${prefix}.vcf.gz.tbi
    echo "" | gzip > ${prefix}.complex_CNV_intervals_to_test.bed.gz
    """
}
