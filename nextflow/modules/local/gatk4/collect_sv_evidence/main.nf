//
// Per-sample PE (discordant pairs), SR (split reads), and SD (site
// depth/BAF, at known SNP sites) evidence collection from a BAM/CRAM.
// Wraps GATK's CollectSVEvidence walker (marked **BETA - WORK IN
// PROGRESS** by GATK itself) -- no nf-core module exists for it. Same
// container as gatk4/svcluster/grouped_sv_cluster -- verified directly
// (`gatk CollectSVEvidence --help`) that the walker and every flag name
// used below exist in that image, so no need for GATK-SV's own separate
// gatk_docker image.
//
// Reproduces wdl/CollectSVEvidence.wdl's RunCollectSVEvidence task.
// site_depth_min_mapq/site_depth_min_baseq are left at GATK-SV's own
// defaults (6/10) via task.ext.args, same pattern as every other tunable
// flag in this pipeline.
//
// sd_locs_vcf (GATK-SV's own resources_hg38.json: sd_locs_vcf ->
// Homo_sapiens_assembly38.dbsnp138.vcf) is a plain, uncompressed VCF with
// a GATK-style .idx sibling index -- not bgzipped/.tbi like every other
// VCF this pipeline handles. -F/--site-depth-locs-vcf takes only the VCF
// path (no separate index flag, confirmed via `gatk CollectSVEvidence
// --help`); the .idx just needs to be staged alongside it for GATK's own
// auto-discovery to find, hence the plain `path` (not `tuple`) input
// below.
//
// Output is NOT indexed here -- `tabix` isn't present in the
// gatk4-main_gcnvkernel container (verified directly: `which tabix` found
// nothing), so indexing is a separate downstream step
// (modules/local/tabix_sv_evidence) in whatever container it's used from,
// same "small chained modules" pattern as
// EXCLUDE_VARIANTS_BY_ID/TABIX elsewhere in this pipeline.
//
process GATK_COLLECT_SV_EVIDENCE {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/b9/b9822b92da68a3e7916072218082e3fa79bebc2f377947c363613adeecd56ec5/data'
        : 'community.wave.seqera.io/library/gatk4-main_gcnvkernel:961440660027ec01'}"

    input:
    tuple val(meta), path(bam), path(bai)
    path(fasta)
    path(fasta_fai)
    path(dict)
    path(sd_locs_vcf)
    path(sd_locs_vcf_idx)

    output:
    tuple val(meta), path("*.pe.txt.gz"), emit: pe
    tuple val(meta), path("*.sr.txt.gz"), emit: sr
    tuple val(meta), path("*.sd.txt.gz"), emit: sd

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    def avail_mem = 3072
    if (!task.memory) {
        log.info('[GATK CollectSVEvidence] Available memory not known - defaulting to 3GB. Specify process memory requirements to change this.')
    }
    else {
        avail_mem = (task.memory.mega * 0.8).intValue()
    }
    """
    gatk --java-options "-Xmx${avail_mem}M -XX:-UsePerfData" \\
        CollectSVEvidence \\
        -I ${bam} \\
        --sample-name ${meta.ped_id} \\
        -F ${sd_locs_vcf} \\
        -SR ${prefix}.sr.txt.gz \\
        -PE ${prefix}.pe.txt.gz \\
        -SD ${prefix}.sd.txt.gz \\
        -R ${fasta} \\
        --read-filter NonZeroReferenceLengthAlignmentReadFilter \\
        ${args}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo "" | gzip > ${prefix}.pe.txt.gz
    echo "" | gzip > ${prefix}.sr.txt.gz
    echo "" | gzip > ${prefix}.sd.txt.gz
    """
}
