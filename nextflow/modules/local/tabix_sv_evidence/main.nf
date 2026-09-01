//
// Index one of GATK's SV evidence files (PE/SR/SD -- see
// modules/local/gatk4/collect_sv_evidence) with tabix. These are 0-based,
// tab-delimited, sorted-by-position files, not VCFs, so they need the
// same explicit `-0 -s1 -b2 -e2` flags GATK-SV's own WDL task uses
// (wdl/CollectSVEvidence.wdl) -- plain `tabix <file>` (see
// modules/local/tabix, used for VCFs) can't auto-detect this format.
// Separate module rather than parameterizing modules/local/tabix, to keep
// that one's signature simple for its one real use case (VCFs).
//
// GATK's CollectSVEvidence walker writes its PE/SR/SD output already
// bgzipped (an htsjdk feature-file writer, same mechanism GATK's VCF
// writers use for .vcf.gz) -- confirmed by wdl/CollectSVEvidence.wdl
// itself, which calls `tabix` directly on the walker's output with no
// separate `bgzip` step first. Only indexing is needed here.
//
// A separate module (not folded into GATK_COLLECT_SV_EVIDENCE's own
// process) because `tabix` isn't present in the gatk4-main_gcnvkernel
// container GATK_COLLECT_SV_EVIDENCE uses (verified directly: `which
// tabix` found nothing in that image) -- a Nextflow process has exactly
// one container, so indexing has to be a separate process using the
// bcftools/htslib container that already provides it (same one
// modules/local/tabix and modules/local/exclude_variants_by_id use).
//
process TABIX_SV_EVIDENCE {
    tag "$meta.id"
    label 'process_single'

    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/0b/0b4d52ca9a56d07be3f78a12af654e5116f5112908dba277e6796fd9dfb83fe5/data' :
        'community.wave.seqera.io/library/bcftools_htslib:1.23.1--9f08ec665533d64a' }"

    input:
    tuple val(meta), path(evidence_file)

    output:
    tuple val(meta), path(evidence_file), path("*.tbi"), emit: evidence

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    tabix -f -0 -s1 -b2 -e2 ${evidence_file}
    """

    stub:
    """
    touch ${evidence_file}.tbi
    """
}
