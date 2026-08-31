//
// Index an already-bgzipped VCF with tabix. A small local module rather
// than nf-core's tabix/tabix (deprecated, hard-fails -- points at
// htslib/bgziptabix instead) or htslib/bgziptabix itself (a general
// compress-and/or-index module with file-type-sniffing branches for
// inputs that may or may not already be compressed -- more machinery than
// a single `tabix <file}` shell-out needs when the input is always
// already bgzipped, as it is here). Same container as EXCLUDE_VARIANTS_BY_ID,
// which already depends on bcftools/htslib.
//
process TABIX {
    tag "$meta.id"
    label 'process_single'

    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/0b/0b4d52ca9a56d07be3f78a12af654e5116f5112908dba277e6796fd9dfb83fe5/data' :
        'community.wave.seqera.io/library/bcftools_htslib:1.23.1--9f08ec665533d64a' }"

    input:
    tuple val(meta), path(vcf)

    output:
    tuple val(meta), path(vcf), path("*.tbi"), emit: vcf

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    tabix ${vcf}
    """

    stub:
    """
    touch ${vcf}.tbi
    """
}
