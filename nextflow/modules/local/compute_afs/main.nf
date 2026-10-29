//
// Sex-aware AN/AC/AF (plus per-sex AN_MALE/AF_FEMALE etc., genotype
// counts and frequencies) with upstream's compute_AFs.py, as
// wdl/ShardedAnnotateVcf.wdl's ComputeAFs does. Outside the PARs on chrX
// and chrY, males count as hemizygous: one allele each, and AC counts
// carrier genotypes. GATK-SV keeps male genotypes in diploid notation
// there (CleanVcf's cleanvcf_preprocess.py), so bcftools +fill-tags would
// count males as diploid.
//
// Unsharded (upstream shards by records), no population file.
//
process COMPUTE_AFS {
    tag "$meta.id"
    label 'process_single'

    container 'us.gcr.io/broad-dsde-methods/gatk-sv/sv-pipeline:2026-06-02-v1.1-24f0fd08'

    input:
    tuple val(meta), path(vcf), path(vcf_index)
    path(ped_file)      // PED_SEXES' ID+sex table: the script reads only ID and sex
    path(par_bed)
    val(chr_x)
    val(chr_y)

    output:
    tuple val(meta), path("${task.ext.prefix ?: meta.id}.vcf.gz"), path("${task.ext.prefix ?: meta.id}.vcf.gz.tbi"), emit: vcf

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    if ("${prefix}.vcf.gz" == vcf.name) {
        error "COMPUTE_AFS: output ${prefix}.vcf.gz would overwrite its input; give it a different meta.id"
    }
    """
    # The script takes the first tab-separated column of each line.
    printf '%s\\t0\\n' ${chr_x} ${chr_y} > allosomes.list

    python /opt/sv-pipeline/05_annotation/scripts/compute_AFs.py ${vcf} stdout \\
        -f ${ped_file} \\
        --par ${par_bed} \\
        --allosomes-list allosomes.list \\
        | bgzip -c > ${prefix}.vcf.gz
    tabix -p vcf ${prefix}.vcf.gz
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo "" | gzip > ${prefix}.vcf.gz
    touch ${prefix}.vcf.gz.tbi
    """
}
