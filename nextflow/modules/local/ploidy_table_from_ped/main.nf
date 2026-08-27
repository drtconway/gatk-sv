//
// Build a per-sample, per-contig ploidy table from a PED file, for
// GATK's SVCluster (and later, genotyping). Wraps GATK-SV's own
// ploidy_table_from_ped.py (vendored in ../../../bin/, no third-party
// dependencies beyond the stdlib) -- see wdl/TasksClusterBatch.wdl's
// CreatePloidyTableFromPed task, which this reproduces.
//
process PLOIDY_TABLE_FROM_PED {
    tag "$meta.id"
    label 'process_single'

    // No third-party deps needed (stdlib only), but reuses the same image
    // as the other vendored-script modules (nextflow/dockerfiles/sv-scripts)
    // rather than pulling a second, otherwise-unneeded image on HPC.
    container 'docker://drtomc/gatk-sv-nf-sv-scripts:0.1.0'

    input:
    tuple val(meta), path(ped)
    path(contig_list)

    output:
    tuple val(meta), path("*.ploidy.tsv"), emit: ploidy_table

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    ploidy_table_from_ped.py \\
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
