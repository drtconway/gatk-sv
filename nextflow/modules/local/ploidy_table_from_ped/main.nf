//
// Build a per-sample, per-contig ploidy table from a PED file, for
// GATK's SVCluster (and later, genotyping). Wraps GATK-SV's own
// ploidy_table_from_ped.py (vendored in ../../../bin/, no third-party
// dependencies beyond the stdlib) -- see wdl/TasksClusterBatch.wdl's
// CreatePloidyTableFromPed task, which this reproduces.
//
// The script is passed in as an explicit `path` input rather than relying
// on Nextflow's bin/-auto-PATH mechanism: that mechanism only adds bin/
// to PATH when it's a sibling of the *invoked* script (baseDir =
// scriptFile.main.parent, unconditionally -- not the launch/current
// directory), so it silently doesn't apply to any workflows/*.nf entry
// point, since bin/ lives in nextflow/, not nextflow/workflows/. This bit
// us for real (exit 127, "command not found") on an HPC run even though
// local -stub-run never caught it -- the stub: block never calls the
// script by name, only script: does.
//
process PLOIDY_TABLE_FROM_PED {
    tag "$meta.id"
    label 'process_single'

    // No third-party deps needed (stdlib only), but reuses the same image
    // as the other vendored-script modules (nextflow/dockerfiles/sv-scripts)
    // rather than pulling a second, otherwise-unneeded image on HPC.
    // Plain image reference, not a docker:// URI: Apptainer/Singularity
    // accept "docker://<image>" specifically to mean "pull from Docker
    // Hub", but plain `docker` itself does not -- it wants the bare
    // reference. A hardcoded "docker://" prefix here broke `-profile
    // docker` outright ("docker: invalid reference format"), only caught
    // by a real (non-stub) local run since -stub-run never invokes the
    // container engine at all.
    container 'drtomc/gatk-sv-nf-sv-scripts:0.1.0'

    input:
    tuple val(meta), path(ped)
    path(contig_list)
    path(script)

    output:
    tuple val(meta), path("*.ploidy.tsv"), emit: ploidy_table

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    python3 ${script} \\
        --ped ${ped} \\
        --contigs ${contig_list} \\
        --out ${prefix}.ploidy.tsv \\
        ${args}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.ploidy.tsv
    """
}
