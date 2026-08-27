//
// Filter a VCF to drop records whose ID appears in an exclusion list,
// then index. The `ID!=@file` bcftools expression needs the exclusion
// file's staged path known at script-render time, which no generic
// nf-core bcftools module input (regions/targets/samples-file) maps onto
// -- those are different bcftools flags entirely -- so this is a small
// local module rather than a repurposed nf-core one.
//
// Used twice in vcfs_cluster_svcluster: pre-clustering (GATK-SV's
// PreparePESRVcfs, wdl/PESRClustering.wdl) and post-clustering (GATK-SV's
// ExcludeIntervalsByEndpoints, wdl/TasksClusterBatch.wdl) -- both are the
// same "bcftools view -i 'ID!=@list' | tabix" shape, just applied at
// different pipeline stages with different exclusion lists.
//
process EXCLUDE_VARIANTS_BY_ID {
    tag "$meta.id"
    label 'process_single'

    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/0b/0b4d52ca9a56d07be3f78a12af654e5116f5112908dba277e6796fd9dfb83fe5/data' :
        'community.wave.seqera.io/library/bcftools_htslib:1.23.1--9f08ec665533d64a' }"

    input:
    tuple val(meta), path(vcf), path(excluded_ids)

    output:
    tuple val(meta), path("*.filtered.vcf.gz"), path("*.filtered.vcf.gz.tbi"), emit: vcf

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    bcftools view -i 'ID!=@${excluded_ids}' ${vcf} -Oz -o ${prefix}.filtered.vcf.gz
    tabix ${prefix}.filtered.vcf.gz
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo "" | gzip > ${prefix}.filtered.vcf.gz
    touch ${prefix}.filtered.vcf.gz.tbi
    """
}
