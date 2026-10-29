//
// Merge a contig's inversion-only and all-variant resolved VCFs, then give
// every record a consistent ID (<cohort>_<SVTYPE>_<contig>_<n>). Reproduces
// wdl/ResolveComplexVariants.wdl's IntegrateResolvedVcfs + RenameVariants.
// rename.py also drops records under 50bp and turns unresolved INV into
// BND.
//
process INTEGRATE_RESOLVED_VCFS {
    tag "${meta.id}"
    label 'process_single'

    container 'us.gcr.io/broad-dsde-methods/gatk-sv/sv-pipeline:2026-06-02-v1.1-24f0fd08'

    input:
    tuple val(meta), path(inv_vcf), path(inv_vcf_index), path(all_vcf), path(all_vcf_index)  // meta.contig, meta.cohort

    output:
    tuple val(meta), path("${task.ext.prefix ?: meta.id}.vcf.gz"), path("${task.ext.prefix ?: meta.id}.vcf.gz.tbi"), emit: vcf

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p tmp
    python /opt/sv-pipeline/04_variant_resolution/scripts/integrate_resolved_vcfs.py \\
        --all-vcf ${all_vcf} \\
        --inv-only-vcf ${inv_vcf} \\
        | bcftools sort --temp-dir ./tmp -Oz -o integrated.vcf.gz
    tabix integrated.vcf.gz

    python /opt/sv-pipeline/04_variant_resolution/scripts/rename.py \\
        --chrom ${meta.contig} --prefix ${meta.cohort} integrated.vcf.gz - \\
        | bgzip > ${prefix}.vcf.gz
    tabix ${prefix}.vcf.gz
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo "" | gzip > ${prefix}.vcf.gz
    touch ${prefix}.vcf.gz.tbi
    """
}
