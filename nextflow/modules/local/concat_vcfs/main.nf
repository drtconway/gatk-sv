//
// Concatenate multiple VCFs into one. Wraps `bcftools concat` -- no
// nf-core module matches this exact shape (nf-core's bcftools/concat
// assumes coordinate-sorted, non-overlapping inputs by default; this
// pipeline's two call sites need --allow-overlaps and --naive
// respectively, via task.ext.args, same as every other per-call-site
// flag difference in this pipeline -- see conf/modules.config). Same
// bcftools_htslib container as TABIX/BGZIP/SCATTER_VCF elsewhere in this
// pipeline. Reproduces wdl/TasksMakeCohortVcf.wdl's ConcatVcfs task,
// minus its sites_only/sort_vcf_list options -- no caller here needs
// them yet (GenerateBatchMetrics.wdl's own two call sites use neither).
//
process CONCAT_VCFS {
    tag "${meta.id}"
    label 'process_single'

    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/0b/0b4d52ca9a56d07be3f78a12af654e5116f5112908dba277e6796fd9dfb83fe5/data' :
        'community.wave.seqera.io/library/bcftools_htslib:1.23.1--9f08ec665533d64a' }"

    input:
    tuple val(meta), path(vcfs), path(vcf_indices)

    output:
    tuple val(meta), path("*.vcf.gz"), path("*.vcf.gz.tbi"), emit: vcf

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    printf '%s\\n' ${vcfs} > vcfs.list

    bcftools concat --no-version ${args} -Oz --file-list vcfs.list > ${prefix}.vcf.gz
    tabix ${prefix}.vcf.gz
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo "" | gzip > ${prefix}.vcf.gz
    touch ${prefix}.vcf.gz.tbi
    """
}
