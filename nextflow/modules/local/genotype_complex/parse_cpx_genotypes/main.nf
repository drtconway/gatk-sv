//
// Merge one contig's RdTest genotypes and apply them to its complex SVs:
// confirm, reclassify or drop each CPX/unresolved-INV record by depth.
// Reproduces wdl/GenotypeCpxCnvs.wdl's MergeMeltedGts + ParseGenotypes.
// The PED is PED_SEXES' ID+sex table (the script reads only ID and sex).
//
process PARSE_CPX_GENOTYPES {
    tag "${meta.id}"
    label 'process_single'

    container 'us.gcr.io/broad-dsde-methods/gatk-sv/sv-pipeline:2026-06-02-v1.1-24f0fd08'

    input:
    tuple val(meta), path(vcf), path(vcf_index), path(intervals), path(genotypes, stageAs: 'genos/*')  // meta.contig
    path(ped_file)

    output:
    tuple val(meta), path("${task.ext.prefix ?: meta.id}.vcf.gz"), path("${task.ext.prefix ?: meta.id}.vcf.gz.tbi"), emit: vcf
    tuple val(meta), path("*.CPXregenotyping_reclassification_table.*.txt.gz"), emit: reclassification_table

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def args = task.ext.args ?: ''  // e.g. --min-ddup-thresh
    """
    mkdir -p tmp
    zcat genos/* | sort -T ./tmp -Vk1,1 -k2,2n -k3,3n -k4,4V -k5,5V | bgzip -c > merged_rd_genos.bed.gz

    python /opt/sv-pipeline/04_variant_resolution/scripts/process_posthoc_cpx_depth_regenotyping.py \\
        --vcf ${vcf} \\
        --intervals ${intervals} \\
        --genotypes merged_rd_genos.bed.gz \\
        ${args} \\
        --ped ${ped_file} \\
        --out out.vcf.gz \\
        --reclassification-table ${prefix}.CPXregenotyping_reclassification_table.${meta.contig}.txt
    gzip ${prefix}.CPXregenotyping_reclassification_table.${meta.contig}.txt

    # Coordinates may have changed
    bcftools sort --temp-dir ./tmp -Oz -o ${prefix}.vcf.gz out.vcf.gz
    tabix ${prefix}.vcf.gz
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo "" | gzip > ${prefix}.vcf.gz
    touch ${prefix}.vcf.gz.tbi
    echo "" | gzip > ${prefix}.CPXregenotyping_reclassification_table.${meta.contig}.txt.gz
    """
}
