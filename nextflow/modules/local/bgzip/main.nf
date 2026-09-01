//
// Compress a plain-text file with bgzip. A small local module rather
// than nf-core's htslib/bgziptabix (a general compress-and/or-index
// module with file-type-sniffing branches this single-purpose call
// doesn't need) -- same reasoning as modules/local/tabix. Same container
// (bcftools/htslib), used here for `bgzip` rather than `tabix`/`bcftools`.
//
process BGZIP {
    tag "$meta.id"
    label 'process_single'

    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/0b/0b4d52ca9a56d07be3f78a12af654e5116f5112908dba277e6796fd9dfb83fe5/data' :
        'community.wave.seqera.io/library/bcftools_htslib:1.23.1--9f08ec665533d64a' }"

    input:
    tuple val(meta), path(infile)

    output:
    tuple val(meta), path("${infile}.gz"), emit: gz

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    bgzip ${infile}
    """

    stub:
    """
    echo "" | gzip > ${infile}.gz
    """
}
