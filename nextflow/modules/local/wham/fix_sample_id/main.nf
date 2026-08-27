//
// Rewrite a Wham VCF's sample column and TAGS INFO field to the pipeline's
// sample_id.
//
// Wham defaults to the BAM's own SM tag for both the VCF sample column and
// the TAGS INFO field (which it also uses to store the sample identifier).
// GATK-SV's own Wham task does this same rewrite inline
// (wdl/Whamg.wdl:162-169) because its svtk standardize_vcf step reads TAGS
// specifically for WHAM VCFs -- a mismatch here means WHAM records get
// attributed to the wrong sample downstream. No nf-core module covers the
// TAGS half (it's Wham-specific), so this is a small local process rather
// than a generic bcftools wrapper.
//
process WHAM_FIX_SAMPLE_ID {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/0b/0b4d52ca9a56d07be3f78a12af654e5116f5112908dba277e6796fd9dfb83fe5/data' :
        'community.wave.seqera.io/library/bcftools_htslib:1.23.1--9f08ec665533d64a' }"

    input:
    tuple val(meta), path(vcf), path(tbi)

    output:
    tuple val(meta), path("*.vcf.gz")     , emit: vcf
    tuple val(meta), path("*.vcf.gz.tbi") , emit: tbi

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo "${meta.id}" > samples.txt

    bcftools reheader --samples samples.txt "${vcf}" \\
        | bcftools view \\
        | sed -e 's/;TAGS=[^;]*;/;TAGS=${meta.id};/' \\
        | bgzip -c > "${prefix}.wham.vcf.gz"

    tabix -p vcf "${prefix}.wham.vcf.gz"
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo "" | gzip > ${prefix}.wham.vcf.gz
    touch ${prefix}.wham.vcf.gz.tbi
    """
}
