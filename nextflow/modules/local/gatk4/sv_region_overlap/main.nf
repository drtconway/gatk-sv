//
// Annotate a VCF's variants with overlap fraction/endpoint-count metrics
// against one or more interval tracks (segdup/repeat-mask), for
// AggregateSVEvidence/AggregateTests downstream. Wraps GATK's
// SVRegionOverlap walker -- no nf-core module exists for it. Reproduces
// wdl/GenerateBatchMetrics.wdl's SVRegionOverlap task.
//
// Different container from every other gatk4/* module here:
// SVRegionOverlap (like AggregateSVEvidence/AggregateDepthEvidence,
// elsewhere in this directory, and TrainSVGenotyping/GenotypeSVs still to
// come) only exists in GATK-SV's own patched GATK build, not the
// bioconda gatk4-main package the gatk4-main_gcnvkernel image (used by
// GATK4_SVCLUSTER, GATK_PRINT_SV_EVIDENCE, etc.) is built from --
// confirmed directly: `gatk --list` on gatk4-main_gcnvkernel has no
// SVRegionOverlap at all, while GATK-SV's own gatk_docker
// (inputs/values/dockers.json) does, at a different GATK version
// (4.6.2.0-92-g6379d28-SNAPSHOT vs bioconda's 4.7.0.0). Using GATK-SV's
// own image directly here rather than waiting on/building a smaller
// substitute.
//
// Plain image reference (not a docker:// URI, not engine-conditional):
// same reasoning as modules/local/ploidy_table_from_ped's own note --
// Apptainer/Singularity accept a bare "registry/path:tag" reference (they
// default to Docker-Hub-style pulling for any unscoped reference, not
// just Docker Hub itself), and plain `docker` requires the bare form
// anyway, so one reference works for both engines here.
//
process GATK_SVREGIONOVERLAP {
    tag "${meta.id}"
    label 'process_single'

    container 'us.gcr.io/broad-dsde-methods/gatk-sv/gatk:mw-gatk-sv-53d5c2d'

    input:
    tuple val(meta), path(vcf), path(vcf_index)
    path(reference_dict)
    tuple path(segdups), path(segdups_index)
    tuple path(rmsk), path(rmsk_index)

    output:
    tuple val(meta), path("*.vcf.gz"), path("*.vcf.gz.tbi"), emit: vcf

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    def avail_mem = 3072
    if (!task.memory) {
        log.info('[GATK SVRegionOverlap] Available memory not known - defaulting to 3GB. Specify process memory requirements to change this.')
    }
    else {
        avail_mem = (task.memory.mega * 0.8).intValue()
    }
    """
    gatk --java-options "-Xmx${avail_mem}M -XX:-UsePerfData" \\
        SVRegionOverlap \\
        -V ${vcf} \\
        -O ${prefix}.vcf.gz \\
        --sequence-dictionary ${reference_dict} \\
        --track-intervals ${segdups} --track-name SEGDUP \\
        --track-intervals ${rmsk} --track-name RMSK \\
        ${args}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo "" | gzip > ${prefix}.vcf.gz
    touch ${prefix}.vcf.gz.tbi
    """
}
