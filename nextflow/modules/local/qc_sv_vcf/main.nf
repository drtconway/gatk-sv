//
// QC of one SV VCF with scripts/qc_genotyped_vcf.sh: site counts, sizes
// and frequencies by SVTYPE, per-sample counts with upstream's IQR outlier
// rule, and trio Mendelian consistency, de novo and transmission rates,
// pooled by stratum, with a GQ sweep. Not in GATK-SV.
//
// meta.filters (optional, e.g. 'PASS,.') restricts the QC to records with
// those FILTER values first. Takes the full PED: trios need parentage, so
// a pedigree edit reruns the QC, which is what it should do.
//
process QC_SV_VCF {
    tag "$meta.id"
    label 'process_single'

    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/0b/0b4d52ca9a56d07be3f78a12af654e5116f5112908dba277e6796fd9dfb83fe5/data' :
        'community.wave.seqera.io/library/bcftools_htslib:1.23.1--9f08ec665533d64a' }"

    input:
    tuple val(meta), path(vcf), path(vcf_index)
    path(ped)
    path(qc_script)   // scripts/qc_genotyped_vcf.sh

    output:
    tuple val(meta), path("${task.ext.prefix ?: meta.id}"), emit: qc

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    // Named after the QC, so summary.txt's first line says which subset it is.
    def input_vcf = meta.filters ? "${prefix}.vcf.gz" : "${vcf}"
    def subset = meta.filters ? "bcftools view --no-version -f '${meta.filters}' ${vcf} -Oz -o ${prefix}.vcf.gz" : ''
    """
    ${subset}
    bash ${qc_script} ${input_vcf} ${ped} ${prefix} > /dev/null
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p ${prefix}
    touch ${prefix}/summary.txt
    """
}
