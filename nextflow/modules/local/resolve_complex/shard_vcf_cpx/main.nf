//
// Group candidate complex-SV records into clusters (svtk vcfcluster,
// sites only, no merging), then split the variant IDs into shards of about
// records_per_shard without splitting a cluster. Reproduces
// wdl/ResolveCpxSv.wdl's ShardVcfCpx plus TasksMakeCohortVcf.wdl's
// ShardVidsForClustering, in one task. No shard lists are emitted for a
// VCF with no records.
//
process SHARD_VCF_CPX {
    tag "${meta.id}"
    label 'process_single'

    container 'us.gcr.io/broad-dsde-methods/gatk-sv/sv-pipeline:2026-06-02-v1.1-24f0fd08'

    input:
    tuple val(meta), path(vcf), path(vcf_index)
    val(dist)               // svtk vcfcluster -d
    val(frac)               // svtk vcfcluster -f
    val(records_per_shard)

    output:
    tuple val(meta), path("*.vids.shard_*.list"), optional: true, emit: vids

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    bcftools view --no-version -G ${vcf} -Oz -o sites_only.vcf.gz
    svtk vcfcluster <(echo "sites_only.vcf.gz") ${prefix}.shard_cpx.vcf \\
        -d ${dist} \\
        -f ${frac} \\
        --single-end \\
        -p candidate_complex_clusters \\
        --svtypes DEL,DUP,INS,INV,BND \\
        --ignore-svtypes \\
        -o 0 \\
        --preserve-header \\
        --preserve-ids \\
        --skip-merge

    # A cluster is never split across shards; a new shard starts once the
    # current one has reached records_per_shard.
    python3 - ${prefix}.shard_cpx.vcf ${prefix} ${records_per_shard} <<'PYEOF'
    import sys
    import pysam

    path, prefix, per_shard = sys.argv[1], sys.argv[2], int(sys.argv[3])
    state = {'shard': -1, 'size': per_shard, 'fout': None}

    def write_cluster(vids):
        if not vids:
            return
        if state['size'] >= per_shard:
            if state['fout'] is not None:
                state['fout'].close()
            state['shard'] += 1
            state['size'] = 0
            state['fout'] = open('%s.vids.shard_%d.list' % (prefix, state['shard']), 'w')
        for vid in vids:
            state['fout'].write(vid + '\\n')
        state['size'] += len(vids)

    current, vids = None, []
    for record in pysam.VariantFile(path):
        cluster = record.info['CLUSTER']
        if cluster != current:
            write_cluster(vids)
            current, vids = cluster, []
        vids.append(record.id)
    write_cluster(vids)
    if state['fout'] is not None:
        state['fout'].close()
    PYEOF
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo stub > ${prefix}.vids.shard_0.list
    """
}
