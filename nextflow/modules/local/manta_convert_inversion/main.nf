//
// Convert Manta's inversion BND pairs (both mates ++, or both --) into
// single INV records flagged INV3/INV5, with Manta's own
// convertInversion.py, as wdl/Manta.wdl does before standardization.
// svtk standardize then gives them SVTYPE=INV with STRANDS ++/--, which
// ResolveComplexVariants' inversion pass selects (SVTYPE="INV"). Without
// it they stay BND and no inversion is ever resolved (README: "Inversions
// were never resolved").
//
// Runs in upstream's Manta image: the script needs python2 and samtools
// (it fetches REF bases with samtools faidx), and the nf-core Manta image
// has no samtools.
//
process MANTA_CONVERT_INVERSION {
    tag "$meta.id"
    label 'process_single'

    container 'us.gcr.io/broad-dsde-methods/gatk-sv/manta:2023-09-14-v0.28.3-beta-3f22f94d'

    input:
    tuple val(meta), path(vcf)
    tuple path(fasta), path(fasta_fai)

    output:
    tuple val(meta), path("${task.ext.prefix ?: meta.id}.manta.inv.vcf.gz"), emit: vcf

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    python2 /usr/local/bin/manta/libexec/convertInversion.py \\
        /usr/local/bin/samtools \\
        ${fasta} \\
        ${vcf} \\
        | bgzip -c > ${prefix}.manta.inv.vcf.gz
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo "" | gzip > ${prefix}.manta.inv.vcf.gz
    """
}
