//
// Per-contig inputs for ResolveComplexVariants' two resolution passes:
// the contig's inversions only (SubsetInversions), and the contig's full
// VCF after the same-breakpoint overlap filter (BreakpointOverlap, which
// drops the weaker of records sharing a breakpoint). Reproduces those two
// wdl/ResolveComplexVariants.wdl tasks, plus the per-contig split that
// upstream gets from CombineBatches' contig scatter.
//
// The SR lists are genome-wide here; the filter reads IDs from their last
// column, so the extra IDs are harmless.
//
process RESOLVE_PREP_CONTIG {
    tag "${meta.id}"
    label 'process_single'

    container 'us.gcr.io/broad-dsde-methods/gatk-sv/sv-pipeline:2026-06-02-v1.1-24f0fd08'

    input:
    tuple val(meta), path(vcf), path(vcf_index)  // meta.contig: the contig to take
    path(bothside_pass_list)
    path(background_fail_list)

    output:
    tuple val(meta), path("${task.ext.prefix ?: meta.id}.inv_only.vcf.gz"), path("${task.ext.prefix ?: meta.id}.inv_only.vcf.gz.tbi"), emit: inv_vcf
    tuple val(meta), path("${task.ext.prefix ?: meta.id}.breakpoint_overlap.vcf.gz"), path("${task.ext.prefix ?: meta.id}.breakpoint_overlap.vcf.gz.tbi"), emit: all_vcf
    tuple val(meta), path("${task.ext.prefix ?: meta.id}.breakpoint_overlap.dropped_records.vcf.gz"), emit: dropped_vcf

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    bcftools view --no-version -r ${meta.contig} ${vcf} -Oz -o contig.vcf.gz
    tabix contig.vcf.gz

    bcftools view --no-version --no-update -i 'INFO/SVTYPE="INV"' contig.vcf.gz -Oz -o ${prefix}.inv_only.vcf.gz
    tabix ${prefix}.inv_only.vcf.gz

    python /opt/sv-pipeline/04_variant_resolution/scripts/overlap_breakpoint_filter.py \\
        contig.vcf.gz \\
        ${bothside_pass_list} \\
        ${background_fail_list} \\
        ${prefix}.breakpoint_overlap.dropped_records.vcf.gz \\
        | bgzip > ${prefix}.breakpoint_overlap.vcf.gz
    tabix ${prefix}.breakpoint_overlap.vcf.gz
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    for f in inv_only breakpoint_overlap breakpoint_overlap.dropped_records; do
        echo "" | gzip > ${prefix}.\$f.vcf.gz
    done
    touch ${prefix}.inv_only.vcf.gz.tbi ${prefix}.breakpoint_overlap.vcf.gz.tbi
    """
}
