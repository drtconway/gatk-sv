//
// Split a genotyped VCF into depth-only and PESR-only records, by
// INFO/ALGORITHMS. Wraps a small bcftools filter -- no nf-core module
// matches this exact shape. Reproduces wdl/GenotypeBatch.wdl's
// SeparateDepthPesr task. Same container as TABIX/BGZIP/SCATTER_VCF/
// CONCAT_VCFS elsewhere in this pipeline.
//
process SEPARATE_DEPTH_PESR {
    tag "${meta.id}"
    label 'process_single'

    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/0b/0b4d52ca9a56d07be3f78a12af654e5116f5112908dba277e6796fd9dfb83fe5/data' :
        'community.wave.seqera.io/library/bcftools_htslib:1.23.1--9f08ec665533d64a' }"

    input:
    tuple val(meta), path(vcf), path(vcf_index)

    output:
    tuple val(meta), path("*.depth.vcf.gz"), path("*.depth.vcf.gz.tbi"), emit: depth_vcf
    tuple val(meta), path("*.pesr.vcf.gz"), path("*.pesr.vcf.gz.tbi"), emit: pesr_vcf

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    bcftools view -i 'INFO/ALGORITHMS=="depth"' ${vcf} -Oz -o ${prefix}.depth.vcf.gz
    tabix ${prefix}.depth.vcf.gz
    bcftools view -i 'INFO/ALGORITHMS!="depth"' ${vcf} -Oz -o ${prefix}.pesr.vcf.gz
    tabix ${prefix}.pesr.vcf.gz
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo "" | gzip > ${prefix}.depth.vcf.gz
    touch ${prefix}.depth.vcf.gz.tbi
    echo "" | gzip > ${prefix}.pesr.vcf.gz
    touch ${prefix}.pesr.vcf.gz.tbi
    """
}
