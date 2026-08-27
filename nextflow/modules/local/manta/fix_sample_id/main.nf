//
// Rewrite a Manta VCF's sample column to the pipeline's sample_id.
//
// Manta has no --sample-name flag: it derives the VCF's sample column
// from the BAM's own SM read-group tag, which may not match sample_id.
// GATK-SV's own Manta task does this same rewrite inline
// (wdl/Manta.wdl:152: `bcftools reheader -s <(echo "sample_id")`) --
// downstream steps that key off the sample name (e.g. SVCluster's
// ploidy-table lookup in format_svtk_vcf_for_gatk.py) fail with a
// KeyError otherwise, since the ploidy table is keyed by sample_id, not
// whatever the BAM's SM tag happens to be. Found via exactly that
// failure on a real run -- nf-core's manta/germline module doesn't do
// this rewrite, and nothing needed the VCF's *internal* sample name to
// match sample_id until SVCluster's ploidy lookup did.
//
// GATK-SV's Manta task also runs Manta output through convertInversion.py
// (converts Manta's inversion-signaling BND pairs into proper INV
// records) before this reheader step -- not implemented here yet. See
// README.md known gaps.
//
process MANTA_FIX_SAMPLE_ID {
    tag "$meta.id"
    label 'process_low'

    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/0b/0b4d52ca9a56d07be3f78a12af654e5116f5112908dba277e6796fd9dfb83fe5/data' :
        'community.wave.seqera.io/library/bcftools_htslib:1.23.1--9f08ec665533d64a' }"

    input:
    tuple val(meta), path(vcf)

    output:
    tuple val(meta), path("*.manta.vcf.gz")    , emit: vcf
    tuple val(meta), path("*.manta.vcf.gz.tbi"), emit: tbi

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo "${meta.id}" > samples.txt

    bcftools reheader --samples samples.txt "${vcf}" \\
        | bcftools view -Oz -o "${prefix}.manta.vcf.gz"

    tabix -p vcf "${prefix}.manta.vcf.gz"
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo "" | gzip > ${prefix}.manta.vcf.gz
    touch ${prefix}.manta.vcf.gz.tbi
    """
}
