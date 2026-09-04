//
// Genotype SVs at a sites-only VCF's sites, for one contig, using the
// RD/PE/SR cutoff tables TrainSVGenotyping produced. Wraps GATK's
// GenotypeSVs walker -- no nf-core module exists for it. Reproduces
// wdl/GenotypeBatch.wdl's GenotypeSVs task, minus its own inline
// per-contig PrintSVEvidence calls (this module's caller runs those
// separately, via GATK_PRINT_SV_EVIDENCE_CONTIG, so rd_file/pe_file/
// sr_file here are already contig-subset -- no -L flag needed in this
// module itself since the evidence files are already filtered).
//
// Same GATK-SV-fork-only container as modules/local/gatk4/
// sv_region_overlap -- see that module's own note.
//
process GATK_GENOTYPESVS {
    tag "${meta.id}"
    label 'process_single'

    container 'us.gcr.io/broad-dsde-methods/gatk-sv/gatk:mw-gatk-sv-53d5c2d'

    input:
    tuple val(meta), path(vcf), path(vcf_index)
    val(contig)
    path(median_coverage)
    tuple val(rd_meta), path(rd_file), path(rd_file_index)
    tuple val(pe_meta), path(pe_file), path(pe_file_index)
    tuple val(sr_meta), path(sr_file), path(sr_file_index)
    path(reference_dict)
    tuple val(ploidy_meta), path(ploidy_table)
    tuple path(depth_exclusion_intervals), path(depth_exclusion_intervals_index)
    tuple path(pesr_exclusion_intervals), path(pesr_exclusion_intervals_index)
    tuple val(rd_table_meta), path(rd_table)
    tuple val(pe_table_meta), path(pe_table)
    tuple val(sr_table_meta), path(sr_table)

    output:
    tuple val(meta), path("*.vcf.gz"), path("*.vcf.gz.tbi"), emit: vcf

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    def avail_mem = 3072
    if (!task.memory) {
        log.info('[GATK GenotypeSVs] Available memory not known - defaulting to 3GB. Specify process memory requirements to change this.')
    }
    else {
        avail_mem = (task.memory.mega * 0.8).intValue()
    }
    """
    gatk --java-options "-Xmx${avail_mem}M -XX:-UsePerfData" \\
        GenotypeSVs \\
        -V ${vcf} \\
        -O ${prefix}.vcf.gz \\
        -L ${contig} \\
        --median-coverage ${median_coverage} \\
        --rd-file ${rd_file} \\
        --discordant-pairs-file ${pe_file} \\
        --split-reads-file ${sr_file} \\
        --sequence-dictionary ${reference_dict} \\
        --ploidy-table ${ploidy_table} \\
        --pesr-exclusion-intervals ${pesr_exclusion_intervals} \\
        --depth-exclusion-intervals ${depth_exclusion_intervals} \\
        --rd-table ${rd_table} \\
        --pe-table ${pe_table} \\
        --sr-table ${sr_table} \\
        ${args}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo "" | gzip > ${prefix}.vcf.gz
    touch ${prefix}.vcf.gz.tbi
    """
}
