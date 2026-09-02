//
// Merge one evidence type (PE or SR) across all samples in a run into
// one panel-wide file. Wraps GATK's PrintSVEvidence walker (marked
// **EXPERIMENTAL FEATURE** by GATK itself) -- no nf-core module exists
// for it. Same container as gatk4/svcluster/collect_sv_evidence --
// verified directly (`gatk PrintSVEvidence --help`) that the walker and
// every flag used below exist in that image.
//
// Reproduces wdl/BatchEvidenceMerging.wdl's MergeEvidence task, minus its
// rename_samples/subset_primary_contigs options: our per-sample PE/SR
// files (modules/local/gatk4/collect_sv_evidence) are already written
// with --sample-name set to meta.ped_id at collection time (not
// sample_id), so there's nothing to rename here -- see
// subworkflows/local/vcfs_merge_evidence's own top-of-file note for why
// that sidesteps the WDL's rename step entirely.
//
// -F takes a `.list` file of paths (one per line), same convention
// GATK's --sample-names also uses -- confirmed via the WDL's own
// `write_lines(files)` usage, not just assumed from --help's wording.
//
process GATK_PRINT_SV_EVIDENCE {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/b9/b9822b92da68a3e7916072218082e3fa79bebc2f377947c363613adeecd56ec5/data'
        : 'community.wave.seqera.io/library/gatk4-main_gcnvkernel:961440660027ec01'}"

    input:
    tuple val(meta), path(evidence_files), path(evidence_file_indices)
    path(dict)

    output:
    tuple val(meta), path("*.txt.gz"), emit: evidence

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    def avail_mem = 3072
    if (!task.memory) {
        log.info('[GATK PrintSVEvidence] Available memory not known - defaulting to 3GB. Specify process memory requirements to change this.')
    }
    else {
        avail_mem = (task.memory.mega * 0.8).intValue()
    }
    """
    printf '%s\\n' ${evidence_files} > evidence.list

    gatk --java-options "-Xmx${avail_mem}M -XX:-UsePerfData" \\
        PrintSVEvidence \\
        -F evidence.list \\
        --sequence-dictionary ${dict} \\
        -O ${prefix}.txt.gz \\
        ${args}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo "" | gzip > ${prefix}.txt.gz
    """
}
