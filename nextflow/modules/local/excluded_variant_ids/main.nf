//
// Reduce a bedtools intersect output to a unique, sorted list of variant
// IDs to exclude. Third step of GATK-SV's PreparePESRVcfs interval
// filtering (wdl/PESRClustering.wdl) -- the `cut -f4 | sort | uniq` half.
// Runs in the same container as bedtools/intersect (whose output this
// consumes) since cut/sort/uniq are coreutils present in any base image,
// not worth a fourth container reference for.
//
process EXCLUDED_VARIANT_IDS {
    tag "$meta.id"
    label 'process_single'

    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bedtools:2.31.1--hf5e1c6e_0' :
        'quay.io/biocontainers/bedtools:2.31.1--hf5e1c6e_0' }"

    input:
    tuple val(meta), path(intersect_bed)

    output:
    tuple val(meta), path("*.excluded_vids.list"), emit: excluded_ids

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    cut -f4 ${intersect_bed} | sort | uniq > ${prefix}.excluded_vids.list
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.excluded_vids.list
    """
}
