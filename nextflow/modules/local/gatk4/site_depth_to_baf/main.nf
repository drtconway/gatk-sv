//
// Merge all samples' site-depth (SD) evidence across a run into one
// panel-wide BAF (B-allele frequency) file. Wraps GATK's SiteDepthtoBAF
// walker (marked **EXPERIMENTAL FEATURE** by GATK itself) -- no nf-core
// module exists for it. Same container as gatk4/svcluster/
// collect_sv_evidence -- verified directly (`gatk SiteDepthtoBAF --help`)
// that the walker and every flag used below exist in that image.
//
// Reproduces wdl/BatchEvidenceMerging.wdl's SDtoBAF task -- GATK-SV's
// default evidence-merging path when no separate BAF files exist (only
// SD, which is all bam_collect_evidence produces; see
// subworkflows/local/vcfs_merge_evidence's own top-of-file note).
//
// Deliberately does NOT reproduce the WDL's rename_samples step (an awk
// rewrite of each SD file's sample column before merging): our SD files
// are already written with --sample-name set to meta.ped_id at
// collection time (modules/local/gatk4/collect_sv_evidence), so there's
// nothing to rename. Also note SiteDepthtoBAF's own --sample-names flag
// (confirmed via --help) only matters for *.baf.bci output, not the
// *.baf.txt.gz format used here -- so it wouldn't have helped with
// renaming even if we needed it.
//
process GATK_SITE_DEPTH_TO_BAF {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/b9/b9822b92da68a3e7916072218082e3fa79bebc2f377947c363613adeecd56ec5/data'
        : 'community.wave.seqera.io/library/gatk4-main_gcnvkernel:961440660027ec01'}"

    input:
    tuple val(meta), path(sd_files), path(sd_file_indices)
    path(dict)
    path(sd_locs_vcf)
    path(sd_locs_vcf_idx)

    output:
    tuple val(meta), path("*.baf.txt.gz"), emit: baf

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    def avail_mem = 3072
    if (!task.memory) {
        log.info('[GATK SiteDepthtoBAF] Available memory not known - defaulting to 3GB. Specify process memory requirements to change this.')
    }
    else {
        avail_mem = (task.memory.mega * 0.8).intValue()
    }
    """
    printf '%s\\n' ${sd_files} > sd.list

    gatk --java-options "-Xmx${avail_mem}M -XX:-UsePerfData" \\
        SiteDepthtoBAF \\
        -F sd.list \\
        --baf-sites-vcf ${sd_locs_vcf} \\
        --sequence-dictionary ${dict} \\
        -O ${prefix}.baf.txt.gz \\
        ${args}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo "" | gzip > ${prefix}.baf.txt.gz
    """
}
