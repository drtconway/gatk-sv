//
// Extract a BED of variant start/end/CHR2-END2 breakpoints from a VCF, for
// intersecting against an exclude-intervals BED. First step of GATK-SV's
// PreparePESRVcfs interval-exclusion filtering (wdl/PESRClustering.wdl) --
// reproduces the `bcftools query | awk | sort` portion; the actual
// intersect + filtering happen in bedtools/intersect and bcftools/view.
//
process VCF_ENDS_BED {
    tag "$meta.id"
    label 'process_single'

    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/0b/0b4d52ca9a56d07be3f78a12af654e5116f5112908dba277e6796fd9dfb83fe5/data' :
        'community.wave.seqera.io/library/bcftools_htslib:1.23.1--9f08ec665533d64a' }"

    input:
    tuple val(meta), path(vcf)

    output:
    tuple val(meta), path("*.ends.bed"), emit: bed

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    bcftools query -f '%CHROM\\t%POS\\t%POS\\t%ID\\t%SVTYPE\\n%CHROM\\t%END\\t%END\\t%ID\\t%SVTYPE\\n%CHR2\\t%END2\\t%END2\\t%ID\\t%SVTYPE\\n' ${vcf} \\
        | awk '\$1!="." && \$2!="."' \\
        | sort -k1,1V -k2,2n -k3,3n \\
        > ${prefix}.ends.bed
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.ends.bed
    """
}
