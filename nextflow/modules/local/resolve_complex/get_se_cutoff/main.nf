//
// Minimum PE support for svtk resolve's single-ender rescan: the RF PE
// cutoff (PEQ, a Poisson p-value) converted to a count, floored at 4.
// Reproduces wdl/ResolveCpxSv.wdl's GetSeCutoff for a single batch (the
// first quartile over one batch's value is that value).
//
// Runs in GATK-SV's sv-pipeline image and calls its bundled
// /opt/sv-pipeline scripts, as the WDL does; nothing here is vendored.
//
process GET_SE_CUTOFF {
    tag "${meta.id}"
    label 'process_single'

    container 'us.gcr.io/broad-dsde-methods/gatk-sv/sv-pipeline:2026-06-02-v1.1-24f0fd08'

    input:
    tuple val(meta), path(rf_cutoffs)

    output:
    tuple val(meta), path("se_cutoff.txt"), emit: cutoff

    when:
    task.ext.when == null || task.ext.when

    script:
    // `awk ... exit` rather than upstream's `| head -n1`: head closing the
    // pipe early can SIGPIPE awk under pipefail.
    """
    PEQ=\$(awk -F '\\t' '\$5 == "PEQ" { print \$2; exit }' ${rf_cutoffs})
    python /opt/sv-pipeline/04_variant_resolution/scripts/convert_poisson_p.py "\$PEQ" \\
        | Rscript -e "cat(max(c(4,floor(quantile(as.numeric(scan('stdin',quiet=T)),probs=0.25)))),sep='\\n')" \\
        > se_cutoff.txt
    """

    stub:
    """
    echo 4 > se_cutoff.txt
    """
}
