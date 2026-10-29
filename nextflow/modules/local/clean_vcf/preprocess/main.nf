//
// First CleanVcf step for one contig: take the contig's records, convert
// svtk -> GATK format (ECN from the ploidy table, SR flags from the lists),
// then cleanvcf_preprocess.py (sex-chromosome genotypes, unknown-sex
// samples, SR flag filters). Reproduces wdl/CleanVcfChromosome.wdl's
// FormatVcfToClean (FormatVcfForGatk.wdl) + CleanVcfPreprocess, unsharded
// (upstream: 5,000 records per shard; a contig here has a few thousand).
//
// Uses upstream's own format_svtk_vcf_for_gatk.py from the image, not this
// pipeline's FORMAT_SVTK_VCF_FOR_GATK, which doesn't take the SR lists
// (and changing it would rerun stage 1).
//
process CLEANVCF_PREPROCESS {
    tag "${meta.id}"
    label 'process_single'

    container 'us.gcr.io/broad-dsde-methods/gatk-sv/sv-pipeline:2026-06-02-v1.1-24f0fd08'

    input:
    tuple val(meta), path(vcf), path(vcf_index)  // meta.contig
    path(ploidy_table)
    path(bothside_pass_list)
    path(background_fail_list)
    path(ped_file)          // PED_SEXES' ID+sex table: the script reads only sex
    val(chr_x)
    val(chr_y)
    path(preprocess_script)  // bin/cleanvcf_preprocess.py (newer than the image's, see bin/README.md)

    output:
    tuple val(meta), path("${task.ext.prefix ?: meta.id}.vcf.gz"), path("${task.ext.prefix ?: meta.id}.vcf.gz.tbi"), emit: vcf

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    bcftools view --no-version -r ${meta.contig} ${vcf} -Oz -o contig.vcf.gz
    tabix contig.vcf.gz

    # --- FormatVcfToClean ---
    python /opt/sv-pipeline/scripts/format_svtk_vcf_for_gatk.py \\
        --vcf contig.vcf.gz \\
        --out formatted.vcf.gz \\
        --ploidy-table ${ploidy_table} \\
        --bothside-pass-list ${bothside_pass_list} \\
        --background-fail-list ${background_fail_list}
    tabix formatted.vcf.gz

    # --- CleanVcfPreprocess ---
    bcftools view --header-only formatted.vcf.gz | grep '^##' > header.txt
    printf '%s\\n' \\
        '##FILTER=<ID=UNRESOLVED,Description="Variant is unresolved">' \\
        '##INFO=<ID=HIGH_SR_BACKGROUND,Number=0,Type=Flag,Description="Variant has high number of SR splits in background samples">' \\
        '##INFO=<ID=BOTHSIDES_SUPPORT,Number=0,Type=Flag,Description="Variant has read-level support for both sides of breakpoint">' \\
        '##INFO=<ID=REVISED_EVENT,Number=0,Type=Flag,Description="Variant has been revised due to a copy number mismatch">' \\
        >> header.txt
    bcftools view --header-only formatted.vcf.gz | grep '^#CHROM' >> header.txt
    bcftools view formatted.vcf.gz | bcftools reheader -h header.txt | bgzip -c > processed.reheader.vcf.gz

    python ${preprocess_script} \\
        -V processed.reheader.vcf.gz \\
        -O ${prefix}.vcf.gz \\
        --chrX ${chr_x} \\
        --chrY ${chr_y} \\
        --fail-list ${background_fail_list} \\
        --pass-list ${bothside_pass_list} \\
        --ped-file ${ped_file}
    tabix -p vcf ${prefix}.vcf.gz
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo "" | gzip > ${prefix}.vcf.gz
    touch ${prefix}.vcf.gz.tbi
    """
}
