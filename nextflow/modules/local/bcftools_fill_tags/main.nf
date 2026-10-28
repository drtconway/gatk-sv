//
// Add INFO/AN, AC and AF, computed from each record's GTs (ploidy-aware:
// haploid GTs count one allele). Not part of GATK-SV's GenotypeBatch;
// added so the published VCFs carry cohort allele frequencies. Same
// container as SEPARATE_DEPTH_PESR.
//
// The output is `${prefix}.vcf.gz`, so the caller must give meta.id a
// value that differs from the input's basename: otherwise the output
// writes through the staged input symlink and truncates the upstream
// file (README: "Then `GATK_TRAINSVGENOTYPING` overwrote its own input").
//
process BCFTOOLS_FILL_TAGS {
    tag "${meta.id}"
    label 'process_single'

    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/0b/0b4d52ca9a56d07be3f78a12af654e5116f5112908dba277e6796fd9dfb83fe5/data' :
        'community.wave.seqera.io/library/bcftools_htslib:1.23.1--9f08ec665533d64a' }"

    input:
    tuple val(meta), path(vcf), path(vcf_index)

    output:
    tuple val(meta), path("${task.ext.prefix ?: meta.id}.vcf.gz"), path("${task.ext.prefix ?: meta.id}.vcf.gz.tbi"), emit: vcf

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    if ("${prefix}.vcf.gz" == vcf.name) {
        error "BCFTOOLS_FILL_TAGS: output name ${prefix}.vcf.gz is the same as the input's; give meta.id a distinct value"
    }
    """
    bcftools +fill-tags ${vcf} -Oz -o ${prefix}.vcf.gz --write-index=tbi -- -t AN,AC,AF
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo "" | gzip > ${prefix}.vcf.gz
    touch ${prefix}.vcf.gz.tbi
    """
}
