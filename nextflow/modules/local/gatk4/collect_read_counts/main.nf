//
// Per-sample binned read-depth (RD) coverage collection from a BAM/CRAM,
// over a fixed set of genome-wide intervals. Wraps GATK's
// CollectReadCounts walker -- verified directly (`gatk CollectReadCounts
// --help`) that it exists in the same gatk4-main_gcnvkernel image already
// used by gatk4/svcluster/grouped_sv_cluster/collect_sv_evidence.
//
// Reproduces wdl/CollectCoverage.wdl's CollectCounts task. Output is TSV
// (not GATK-SV's default HDF5) to match GATK-SV's own choice -- their
// downstream evidence-merging step (BatchEvidenceMerging, not yet
// implemented here) expects TSV.
//
// Output is NOT bgzipped here -- `bgzip` isn't present in the
// gatk4-main_gcnvkernel container (verified directly: only `sed` was
// found, not `bgzip`/`tabix`), so compression is a separate downstream
// step (modules/local/bgzip) in a container that has it. Same reasoning
// as modules/local/gatk4/collect_sv_evidence's split with
// modules/local/tabix_sv_evidence.
//
process GATK_COLLECT_READ_COUNTS {
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
    path(preprocessed_intervals)

    output:
    tuple val(meta), path("*.counts.tsv"), emit: counts

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    def avail_mem = 3072
    if (!task.memory) {
        log.info('[GATK CollectReadCounts] Available memory not known - defaulting to 3GB. Specify process memory requirements to change this.')
    }
    else {
        avail_mem = (task.memory.mega * 0.8).intValue()
    }
    """
    gatk --java-options "-Xmx${avail_mem}M -XX:-UsePerfData" \\
        CollectReadCounts \\
        -L ${preprocessed_intervals} \\
        -I ${bam} \\
        -R ${fasta} \\
        --format TSV \\
        --interval-merging-rule OVERLAPPING_ONLY \\
        -O ${prefix}.counts.tsv \\
        ${args}

    sed -ri "s/@RG\\tID:GATKCopyNumber\\tSM:.+/@RG\\tID:GATKCopyNumber\\tSM:${meta.ped_id}/g" ${prefix}.counts.tsv
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.counts.tsv
    """
}
