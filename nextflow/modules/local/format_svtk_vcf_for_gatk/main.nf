//
// Convert a caller's svtk-format VCF to the format GATK's SVCluster
// expects. Wraps GATK-SV's own format_svtk_vcf_for_gatk.py (vendored in
// ../../../bin/) -- see wdl/PESRClustering.wdl's PreparePESRVcfs task,
// which this reproduces the format-conversion half of (interval/size
// filtering is a separate step, see modules/local/filter_pesr_vcf).
//
// Known gap: current pysam, not GATK-SV's pinned pysam==0.15.4 --
// see bin/README.md for the BND/CTX record risk this may carry.
//
// Uses our own image (nextflow/dockerfiles/sv-scripts) with pysam + tabix
// baked in, built and pushed to Docker Hub -- see that Dockerfile for
// build/push instructions and version history.
//
// The script is passed in as an explicit `path` input rather than relying
// on Nextflow's bin/-auto-PATH mechanism -- see the note in
// modules/local/ploidy_table_from_ped/main.nf for why (this bit us for
// real on an HPC run: exit 127, "command not found").
//
process FORMAT_SVTK_VCF_FOR_GATK {
    tag "$meta.id"
    label 'process_single'

    // Plain image reference, not a docker:// URI -- see the note in
    // modules/local/ploidy_table_from_ped/main.nf.
    container 'drtomc/gatk-sv-nf-sv-scripts:0.1.0'

    input:
    tuple val(meta), path(vcf)
    tuple val(meta2), path(ploidy_table)
    path(script)

    output:
    tuple val(meta), path("*.formatted.vcf.gz"), emit: vcf

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    python3 ${script} \\
        --vcf ${vcf} \\
        --out ${prefix}.formatted.vcf.gz \\
        --ploidy-table ${ploidy_table}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo "" | gzip > ${prefix}.formatted.vcf.gz
    """
}
