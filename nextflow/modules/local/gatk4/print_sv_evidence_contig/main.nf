//
// Subset one already-merged, panel-wide evidence file (RD/PE/SR) to a
// single contig. Wraps GATK's PrintSVEvidence walker -- distinct from
// modules/local/gatk4/print_sv_evidence (which merges *multiple
// per-sample* files into one panel-wide file, using -F <list-file>):
// this module takes exactly one already-merged file and filters it with
// -L <contig>, matching wdl/GenotypeBatch.wdl's GenotypeSVs task's own
// per-contig PrintSVEvidence calls exactly (single --evidence-file, not
// a list). Same underlying container/walker as print_sv_evidence -- see
// that module's own note on verification.
//
process GATK_PRINT_SV_EVIDENCE_CONTIG {
    tag "${meta.id}"
    label 'process_single'

    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/b9/b9822b92da68a3e7916072218082e3fa79bebc2f377947c363613adeecd56ec5/data'
        : 'community.wave.seqera.io/library/gatk4-main_gcnvkernel:961440660027ec01'}"

    input:
    tuple val(meta), path(evidence_file), path(evidence_file_index)
    path(dict)
    val(contig)

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
    gatk --java-options "-Xmx${avail_mem}M -XX:-UsePerfData" \\
        PrintSVEvidence \\
        --evidence-file ${evidence_file} \\
        --sequence-dictionary ${dict} \\
        -L ${contig} \\
        -O ${prefix}.txt.gz \\
        ${args}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo "" | gzip > ${prefix}.txt.gz
    """
}
