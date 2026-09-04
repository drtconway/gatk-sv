//
// Extract each depth-genotyped variant's RD_MCR (median coverage ratio)
// FORMAT field into a flat TSV, for RegenotypeCNVs (a downstream
// re-genotyping QC step, not yet built here). Wraps GATK-SV's own
// extract_format_table.py (vendored in ../../../bin/, see
// wdl/GenotypeBatch.wdl's GenerateRegenoCoverageMedians task, which this
// reproduces). Only depends on pysam -- reuses the sv-scripts image
// already used by ploidy_table_from_ped/format_svtk_vcf_for_gatk/
// format_gatk_vcf_for_svtk, rather than pulling a second image.
//
// The script is passed in as an explicit `path` input rather than relying
// on Nextflow's bin/-auto-PATH mechanism -- see
// modules/local/ploidy_table_from_ped/main.nf's own note on why that
// mechanism never fires for any workflows/*.nf entry point here.
//
process GENERATE_REGENO_COVERAGE_MEDIANS {
    tag "${meta.id}"
    label 'process_single'

    container 'drtomc/gatk-sv-nf-sv-scripts:0.1.0'

    input:
    tuple val(meta), path(vcf), path(vcf_index)
    path(script)

    output:
    tuple val(meta), path("*.tsv.gz"), emit: regeno_coverage_medians

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    python3 ${script} \\
        --format-field RD_MCR \\
        --id-column cnvID \\
        --vcf ${vcf} \\
        --out ${prefix}.tsv

    gzip ${prefix}.tsv
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo "" | gzip > ${prefix}.tsv.gz
    """
}
