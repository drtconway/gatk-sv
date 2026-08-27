//
// Derive a bedtools-style genome file (contig, length) from a FASTA
// index, for bedtools intersect's -g flag. Same one-liner GATK-SV's own
// tasks use inline (e.g. wdl/PESRClustering.wdl's PreparePESRVcfs:
// `cut -f1,2 reference_fasta_fai > genome.file`) -- a real process here
// rather than inline Groovy so the output lands in the task work
// directory, not wherever the source .fai happens to live (which may not
// even be writable, and shouldn't be written into regardless).
//
process GENOME_FILE {
    tag "$meta.id"
    label 'process_single'

    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/0b/0b4d52ca9a56d07be3f78a12af654e5116f5112908dba277e6796fd9dfb83fe5/data' :
        'community.wave.seqera.io/library/bcftools_htslib:1.23.1--9f08ec665533d64a' }"

    input:
    tuple val(meta), path(fasta_fai)

    output:
    tuple val(meta), path("genome.file"), emit: genome_file

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    cut -f1,2 ${fasta_fai} > genome.file
    """

    stub:
    """
    touch genome.file
    """
}
