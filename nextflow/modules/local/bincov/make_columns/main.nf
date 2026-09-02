//
// Reshape one sample's CollectReadCounts output to the panel-wide bin
// grid (from BINCOV_SET_BINS), producing one column ready to be pasted
// alongside every other sample's. Reproduces
// wdl/MakeBincovMatrix.wdl's MakeBincovMatrixColumns task.
//
// Same container/tool set as BINCOV_SET_BINS -- see that module's note.
//
process BINCOV_MAKE_COLUMNS {
    tag "${meta.id}"
    label 'process_single'

    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/0b/0b4d52ca9a56d07be3f78a12af654e5116f5112908dba277e6796fd9dfb83fe5/data' :
        'community.wave.seqera.io/library/bcftools_htslib:1.23.1--9f08ec665533d64a' }"

    input:
    tuple val(meta), path(count_file)
    val(binsize)
    path(bin_locs)

    output:
    tuple val(meta), path("*.RD.txt.gz"), emit: bincov_column

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    firstchar=\$(zcat -f ${count_file} | head -c 1)
    if [ "\$firstchar" == '@' ]; then
        shift=1
    else
        shift=0
    fi

    printf '#Chr\\tStart\\tEnd\\t%s\\n' "${meta.ped_id}" > tmp.bed
    zcat -f ${count_file} \\
        | sed '/^@/d' \\
        | sed '/^CONTIG\tSTART\tEND\tCOUNT\$/d' \\
        | sed '/^#/d' \\
        | awk -v x="\${shift}" -v b="${binsize}" 'BEGIN{OFS="\\t"}{\$2=\$2-x; if (\$3-\$2==b) print \$0}' \\
        >> tmp.bed

    if ! cut -f1-3 tmp.bed | cmp <(bgzip -cd ${bin_locs}); then
        echo "${count_file} has different intervals than ${bin_locs}" >&2
        exit 1
    fi
    cut -f4- tmp.bed | bgzip -c > ${prefix}.RD.txt.gz
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo "" | gzip > ${prefix}.RD.txt.gz
    """
}
