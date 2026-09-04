//
// Split a VCF into shards of at most `records_per_shard` records each, for
// per-shard annotation (FormatVcf/SVRegionOverlap/AggregateSVEvidence/
// AggregateDepthEvidence) in vcfs_generate_batch_metrics. Wraps
// `bcftools +scatter` -- no nf-core module exists for this exact
// shape. Reproduces wdl/TasksMakeCohortVcf.wdl's ScatterVcf task,
// including its empty-VCF fallback.
//
// Same bcftools_htslib container as TABIX/BGZIP/EXCLUDE_VARIANTS_BY_ID
// elsewhere in this pipeline -- verified directly (`bcftools +scatter -h`)
// that the plugin exists and every flag used below is real.
//
// The WDL's own placeholder-then-overwrite behavior is deliberate, not
// incidental, and reproduced here: `bcftools +scatter` numbers shards
// starting at 0 and produces ZERO output files for a genuinely
// empty (0-record) input VCF -- confirmed empirically, not assumed from
// the WDL's comment alone. Pre-creating "${prefix}.0.vcf.gz" (a
// header-only shard) before running +scatter means a genuinely empty
// input still yields exactly one (empty) shard rather than zero; for a
// non-empty input, +scatter's own first real shard (also named
// "${prefix}.0...") simply overwrites this placeholder, which is
// harmless -- confirmed directly, not left to chance.
//
process SCATTER_VCF {
    tag "${meta.id}"
    label 'process_single'

    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/0b/0b4d52ca9a56d07be3f78a12af654e5116f5112908dba277e6796fd9dfb83fe5/data' :
        'community.wave.seqera.io/library/bcftools_htslib:1.23.1--9f08ec665533d64a' }"

    input:
    tuple val(meta), path(vcf), path(vcf_index)
    val(records_per_shard)

    output:
    tuple val(meta), path("*.shard_*.vcf.gz"), emit: shards

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def threads = task.cpus ?: 1
    """
    bcftools view -h ${vcf} | bgzip -c > "${prefix}.0.vcf.gz"
    bcftools +scatter ${vcf} -o . -O z -p "${prefix}." --threads ${threads} -n ${records_per_shard}

    ls "${prefix}".*.vcf.gz | sort -k1,1V > vcfs.list
    i=0
    while read -r VCF; do
        shard_no=\$(printf %06d \$i)
        mv "\$VCF" "${prefix}.shard_\${shard_no}.vcf.gz"
        i=\$((i + 1))
    done < vcfs.list
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo "" | gzip > ${prefix}.shard_000000.vcf.gz
    """
}
