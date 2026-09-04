//
// Filter one caller's sites VCF to RF-score-passing (score >= 0.5)
// records, reclassify INV/BND/INS-scoring passers as BND, rewrite
// SR-based breakpoint coordinates, and annotate each surviving record
// with which evidence types (BAF/PE/RD/SR) supported it. Wraps GATK-SV's
// own rewrite_SR_coords.py/annotate_RF_evidence.py (vendored in
// ../../../bin/) plus an inline bcftools/python filtering pipeline.
// Reproduces wdl/FilterBatchSites.wdl's FilterAnnotateVcf task exactly,
// including its own inline BND END2/CHR2 backfill step.
//
// New container (dockerfiles/filter-annotate-vcf): this task's own inline
// shell pipeline needs bcftools/bgzip, and rewrite_SR_coords.py/
// annotate_RF_evidence.py need pysam/numpy/pandas -- neither of this
// pipeline's existing images (bcftools_htslib; aggregate-tests, which has
// the Python packages but no bcftools) covers both, and a single Nextflow
// process has exactly one container, so this combines both rather than
// splitting one WDL task's single inline pipeline across two processes.
//
// The scripts are passed in as explicit `path` inputs rather than
// relying on Nextflow's bin/-auto-PATH mechanism -- see
// modules/local/ploidy_table_from_ped/main.nf's own note on why that
// mechanism never fires for any workflows/*.nf entry point here.
//
process FILTER_ANNOTATE_VCF {
    tag "${meta.id}"
    label 'process_single'

    container 'drtomc/gatk-sv-nf-filter-annotate-vcf:0.1.0'

    input:
    tuple val(meta), path(vcf), path(vcf_index)
    tuple val(metrics_meta), path(metrics)
    tuple val(scores_meta), path(scores)
    tuple val(cutoffs_meta), path(cutoffs)
    path(rewrite_sr_coords_script)
    path(annotate_rf_evidence_script)

    output:
    tuple val(meta), path("*.with_evidence.vcf.gz"), emit: annotated_vcf

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    cat \\
        <(sed -e '1d' ${scores} | fgrep -e DEL -e DUP | awk '(\$3!="NA" && \$3>=0.5)' | cut -f1 | fgrep -w -f - <(zcat ${vcf})) \\
        <(sed -e '1d' ${scores} | fgrep -e INV -e BND -e INS | awk '(\$3!="NA" && \$3>=0.5)' | cut -f1 | fgrep -w -f - <(zcat ${vcf}) | sed -e 's/SVTYPE=DEL/SVTYPE=BND/' -e 's/SVTYPE=DUP/SVTYPE=BND/' -e 's/<DEL>/<BND>/' -e 's/<DUP>/<BND>/') \\
        | cat <(sed -n -e '/^#/p' <(zcat ${vcf})) - \\
        | bcftools sort -Oz -o filtered.vcf.gz

    python3 <<CODE
import pysam

with pysam.VariantFile("filtered.vcf.gz", 'r') as vcf_in, pysam.VariantFile("filtered.updated_bnds.vcf.gz", 'w', header=vcf_in.header) as vcf_out:
    for record in vcf_in:
        if record.info.get('SVTYPE') == 'BND' and 'END2' not in record.info:
            record.info['END2'] = record.stop
            record.stop = record.pos
        if record.info.get('SVTYPE') == 'BND' and 'CHR2' not in record.info:
            record.info['CHR2'] = record.chrom
        vcf_out.write(record)
CODE

    python3 ${rewrite_sr_coords_script} filtered.updated_bnds.vcf.gz ${metrics} ${cutoffs} stdout \\
        | bcftools sort -Oz -o filtered.corrected_coords.vcf.gz

    python3 ${annotate_rf_evidence_script} filtered.corrected_coords.vcf.gz ${scores} ${prefix}.with_evidence.vcf
    bgzip ${prefix}.with_evidence.vcf
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo "" | gzip > ${prefix}.with_evidence.vcf.gz
    """
}
