//
// Context-aware re-clustering of an already-clustered SV VCF, stratifying
// by genomic context tracks (e.g. simple repeats, segmental duplications,
// RepeatMasker). Wraps GATK's GroupedSVCluster walker -- no nf-core module
// exists for it. Same container as gatk4/svcluster (GATK-SV's own WDL
// runs both from the same gatk_docker image, no separate image param for
// GroupedSVCluster -- it's a different walker in the same GATK jar, not a
// separate tool).
//
// Reproduces wdl/CombineBatches.wdl's GroupedSVClusterTask, called twice
// (part 1 with a "coarse" config, part 2 with a "fine" config) by
// subworkflows/local/vcfs_combine_batches -- see that subworkflow for
// where the two rounds' config/stratification files come from
// (clustering_config_part1/2, stratification_config_part1/2, all GATK-SV
// public-bucket resources, same pattern as every other static resource
// so far).
//
process GATK4_GROUPEDSVCLUSTER {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/b9/b9822b92da68a3e7916072218082e3fa79bebc2f377947c363613adeecd56ec5/data'
        : 'community.wave.seqera.io/library/gatk4-main_gcnvkernel:961440660027ec01'}"

    input:
    tuple val(meta), path(vcf), path(vcf_index)
    path(ploidy_table)
    path(fasta)
    path(fasta_fai)
    path(dict)
    path(clustering_config)
    path(stratification_config)
    path(track_bed_files)
    val(track_names)

    output:
    tuple val(meta), path("*.vcf.gz")    , emit: vcf
    tuple val(meta), path("*.vcf.gz.tbi"), emit: vcf_index

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def track_args = track_bed_files.withIndex().collect { bed, i ->
        "--track-intervals ${bed} --track-name ${track_names[i]}"
    }.join(' ')

    def avail_mem = 3072
    if (!task.memory) {
        log.info('[GATK GroupedSVCluster] Available memory not known - defaulting to 3GB. Specify process memory requirements to change this.')
    }
    else {
        avail_mem = (task.memory.mega * 0.8).intValue()
    }
    """
    gatk --java-options "-Xmx${avail_mem}M -XX:-UsePerfData" \\
        GroupedSVCluster \\
        --reference ${fasta} \\
        --ploidy-table ${ploidy_table} \\
        -V ${vcf} \\
        -O ${prefix}.vcf.gz \\
        --clustering-config ${clustering_config} \\
        --stratify-config ${stratification_config} \\
        ${track_args} \\
        --tmp-dir . \\
        ${args}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo "" | gzip > ${prefix}.vcf.gz
    touch ${prefix}.vcf.gz.tbi
    """
}
