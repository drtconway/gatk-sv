//
// Drop VCF records on any contig not in primary_contigs_list (ALT/decoy/
// unplaced contigs, etc). A generic robustness net, not caller-specific:
// GATK-SV restricts *calling* to primary contigs per-caller (e.g. Manta's
// --callRegions, Wham's -c), but a caller can still emit an ALT-contig
// record depending on how it behaves near contig boundaries. Ploidy-table
// lookups downstream (format_svtk_vcf_for_gatk.py) are keyed only by the
// contigs in primary_contigs_list -- a record on any other contig fails
// with a KeyError there (this module exists because a real Manta run on
// real data hit exactly that on an ALT contig: chr1_KI270706v1_random).
//
process VCF_PRIMARY_CONTIGS_ONLY {
    tag "$meta.id"
    label 'process_low'

    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/0b/0b4d52ca9a56d07be3f78a12af654e5116f5112908dba277e6796fd9dfb83fe5/data' :
        'community.wave.seqera.io/library/bcftools_htslib:1.23.1--9f08ec665533d64a' }"

    input:
    tuple val(meta), path(vcf)
    path(contig_list)

    output:
    tuple val(meta), path("*.primary_contigs.vcf.gz")    , emit: vcf
    tuple val(meta), path("*.primary_contigs.vcf.gz.tbi"), emit: tbi

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    tabix -p vcf "${vcf}" 2>/dev/null || true

    CONTIGS=\$(paste -sd, "${contig_list}")
    bcftools view -t "\${CONTIGS}" "${vcf}" -Oz -o "${prefix}.primary_contigs.vcf.gz"
    tabix -p vcf "${prefix}.primary_contigs.vcf.gz"
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo "" | gzip > ${prefix}.primary_contigs.vcf.gz
    touch ${prefix}.primary_contigs.vcf.gz.tbi
    """
}
