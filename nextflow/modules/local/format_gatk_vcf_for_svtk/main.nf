//
// Convert a clustered GATK-format VCF back to svtk format, for downstream
// GATK-SV-compatible tooling. Wraps GATK-SV's own
// format_gatk_vcf_for_svtk.py (vendored in ../../../bin/) -- see
// wdl/TasksClusterBatch.wdl's GatkToSvtkVcf task, which this reproduces.
//
// Known gap: current pysam, not GATK-SV's pinned pysam==0.15.4 -- see
// bin/README.md. Uses our own image (nextflow/dockerfiles/sv-scripts)
// with pysam + tabix baked in -- see that Dockerfile for build/push
// instructions and version history.
//
// The script is passed in as an explicit `path` input rather than relying
// on Nextflow's bin/-auto-PATH mechanism -- see the note in
// modules/local/ploidy_table_from_ped/main.nf for why (this bit us for
// real on an HPC run: exit 127, "command not found").
//
process FORMAT_GATK_VCF_FOR_SVTK {
    tag "$meta.id"
    label 'process_single'

    // Plain image reference, not a docker:// URI -- see the note in
    // modules/local/ploidy_table_from_ped/main.nf.
    container 'drtomc/gatk-sv-nf-sv-scripts:0.1.0'

    input:
    tuple val(meta), path(vcf)
    path(contig_list)
    val(source)
    path(script)

    output:
    tuple val(meta), path("*.svtk.vcf.gz"), path("*.svtk.vcf.gz.tbi"), emit: vcf

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    python3 ${script} \\
        --vcf ${vcf} \\
        --out ${prefix}.svtk.vcf.gz \\
        --source ${source} \\
        --contigs ${contig_list} \\
        ${args}
    tabix ${prefix}.svtk.vcf.gz
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo "" | gzip > ${prefix}.svtk.vcf.gz
    touch ${prefix}.svtk.vcf.gz.tbi
    """
}
