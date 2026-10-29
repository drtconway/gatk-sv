//
// Depth-genotype one chunk of complex-SV CNV intervals with RdTest, on the
// panel-wide read-count matrix subset to the chunk's regions
// (PrintSVEvidence). Reproduces wdl/GenotypeCpxCnvsPerBatch.wdl's
// RdTestGenotype. The PED is PED_SEXES' ID+sex table: RdTest only reads
// the family columns for plotting and de novo mode, both off here.
//
process RDTEST_GENOTYPE {
    tag "${meta.id}"
    label 'process_single'

    container 'us.gcr.io/broad-dsde-methods/gatk-sv/sv-pipeline:2026-06-02-v1.1-24f0fd08'

    input:
    tuple val(meta), path(bed)
    tuple path(coverage_file), path(coverage_file_index)
    path(median_file)
    path(ped_file)
    tuple path(bin_exclude), path(bin_exclude_index)
    path(gt_cutoffs)        // TrainSVGenotyping's rd_table
    path(reference_dict)
    val(n_bins)

    output:
    tuple val(meta), path("${task.ext.prefix ?: meta.id}.rd_genos.bed.gz"), emit: genotypes

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def avail_mem = task.memory ? (task.memory.mega * 0.8).intValue() : 3072
    // The bin_exclude index is staged with it, so upstream's
    // `tabix -p bed bin_exclude` is skipped (it would refuse to overwrite).
    """
    head -1 ${median_file} | sed -e 's/\\t/\\n/g' > samples.list

    # `|| true`: grep exits 1 on an empty chunk, fatal under pipefail
    # (upstream runs this line before `set -o pipefail`).
    (grep -v "^#" ${bed} || true) | sort -T . -k1,1V -k2,2n | bedtools merge -i stdin -d 1000000 > merged.bed
    if [ -s merged.bed ]; then
        java -Xmx${avail_mem}M -jar \$GATK_JAR PrintSVEvidence \\
            --sequence-dictionary ${reference_dict} \\
            --evidence-file ${coverage_file} \\
            -L merged.bed \\
            -O local.RD.txt.gz
        tabix -f -0 -s1 -b2 -e3 local.RD.txt.gz
    else
        touch local.RD.txt
        bgzip local.RD.txt
        tabix -p bed local.RD.txt.gz
    fi

    Rscript /opt/RdTest/RdTest.R \\
        -b ${bed} \\
        -c local.RD.txt.gz \\
        -m ${median_file} \\
        -f ${ped_file} \\
        -n ${prefix} \\
        -w samples.list \\
        -i ${n_bins} \\
        -r ${gt_cutoffs} \\
        -y ${bin_exclude} \\
        -g TRUE

    if [ -f "${prefix}.geno" ] && [ -f "${prefix}.gq" ]; then
        python /opt/sv-pipeline/04_variant_resolution/scripts/merge_RdTest_genotypes.py \\
            ${prefix}.geno ${prefix}.gq rd.geno.cnv.bed
        sort -T . -k1,1V -k2,2n rd.geno.cnv.bed | uniq | bgzip -c > ${prefix}.rd_genos.bed.gz
    else
        # RdTest writes nothing for an empty chunk
        echo "" | bgzip -c > ${prefix}.rd_genos.bed.gz
    fi
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo "" | gzip > ${prefix}.rd_genos.bed.gz
    """
}
