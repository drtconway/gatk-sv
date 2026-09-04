//
// Train genotyping cutoff tables (RD/PE/SR) from a sites-only VCF plus
// panel-wide PE/SR/RD evidence, using rf_cutoffs (from
// FilterBatchSites/AdjudicateSV) to set PE/SR quality thresholds. Wraps
// GATK's TrainSVGenotyping walker -- no nf-core module exists for it.
// Reproduces wdl/GenotypeBatch.wdl's TrainSVGenotyping task.
//
// Same GATK-SV-fork-only container as modules/local/gatk4/
// sv_region_overlap -- see that module's own note on why this differs
// from every other gatk4/* module here (gatk4-main_gcnvkernel doesn't
// have this walker at all; confirmed directly, same as
// AggregateSVEvidence/AggregateDepthEvidence/SVRegionOverlap).
//
process GATK_TRAINSVGENOTYPING {
    tag "${meta.id}"
    label 'process_single'

    container 'us.gcr.io/broad-dsde-methods/gatk-sv/gatk:mw-gatk-sv-53d5c2d'

    input:
    tuple val(meta), path(vcf), path(vcf_index)
    path(training_intervals)
    path(median_coverage)
    tuple val(rd_meta), path(rd_file), path(rd_file_index)
    tuple val(pe_meta), path(pe_file), path(pe_file_index)
    tuple val(sr_meta), path(sr_file), path(sr_file_index)
    path(reference_dict)
    tuple val(ploidy_meta), path(ploidy_table)
    tuple path(depth_exclusion_intervals), path(depth_exclusion_intervals_index)
    tuple path(pesr_exclusion_intervals), path(pesr_exclusion_intervals_index)
    tuple val(cutoffs_meta), path(rf_cutoffs)

    output:
    tuple val(meta), path("*.rd_geno_params.tsv"), emit: rd_table
    tuple val(meta), path("*.pe_geno_params.tsv"), emit: pe_table
    tuple val(meta), path("*.sr_geno_params.tsv"), emit: sr_table

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    def avail_mem = 3072
    if (!task.memory) {
        log.info('[GATK TrainSVGenotyping] Available memory not known - defaulting to 3GB. Specify process memory requirements to change this.')
    }
    else {
        avail_mem = (task.memory.mega * 0.8).intValue()
    }
    """
    PEQ=\$(awk -F '\\t' '{if ( \$5=="PEQ") print \$2 }' ${rf_cutoffs})
    SRQ=\$(awk -F '\\t' '{if ( \$5=="SRQ") print \$2 }' ${rf_cutoffs})

    gatk --java-options "-Xmx${avail_mem}M -XX:-UsePerfData" \\
        TrainSVGenotyping \\
        -XL chrX -XL chrY \\
        -V ${vcf} \\
        --training-intervals ${training_intervals} \\
        -O ${prefix}.vcf.gz \\
        --median-coverage ${median_coverage} \\
        --rd-file ${rd_file} \\
        --split-reads-file ${sr_file} \\
        --discordant-pairs-file ${pe_file} \\
        --sequence-dictionary ${reference_dict} \\
        --ploidy-table ${ploidy_table} \\
        --depth-exclusion-intervals ${depth_exclusion_intervals} \\
        --pesr-exclusion-intervals ${pesr_exclusion_intervals} \\
        --pe-quality \${PEQ} \\
        --sr-quality \${SRQ} \\
        --output-dir ./ \\
        --output-name ${prefix} \\
        ${args}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.rd_geno_params.tsv ${prefix}.pe_geno_params.tsv ${prefix}.sr_geno_params.tsv
    """
}
