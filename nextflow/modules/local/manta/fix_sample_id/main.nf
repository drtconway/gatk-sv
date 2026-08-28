//
// Rewrite a Manta VCF's sample column to the pipeline's ped_id.
//
// Manta has no --sample-name flag: it derives the VCF's sample column
// from the BAM's own SM read-group tag, which may not match anything the
// pipeline uses. GATK-SV's own Manta task does this same rewrite inline
// (wdl/Manta.wdl:152: `bcftools reheader -s <(echo "sample_id")`) --
// downstream steps that key off the sample name (e.g. SVCluster's
// ploidy-table lookup in format_svtk_vcf_for_gatk.py) fail with a
// KeyError otherwise, since the ploidy table is keyed by whatever
// individual_id the pedigree file uses.
//
// Uses ped_id, not sample_id (meta.id): GATK-SV assumes sample_id and the
// PED file's individual_id are the same string, but our real data doesn't
// -- sample_id is often a lab/external reference (e.g.
// "23PS00385-RDN0206"), while the pedigree file uses a different internal
// convention (e.g. "RDN0206-00"). The ploidy table is keyed by whatever
// the PED file's individual_id is, so the VCF's sample column must match
// that, not sample_id -- a real KeyError from exactly this mismatch (using
// sample_id here) is what surfaced the need for ped_id at all. Output
// filenames stay keyed on sample_id/meta.id for readability; only the
// VCF's *internal* sample column uses ped_id.
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
    echo "${meta.ped_id}" > samples.txt

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
