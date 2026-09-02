//
// Paste the bin-location column (BINCOV_SET_BINS) and every sample's
// count column (BINCOV_MAKE_COLUMNS) together, column-wise, into one
// panel-wide RD/bincov matrix. Reproduces wdl/MakeBincovMatrix.wdl's
// ZPaste task -- including its named-pipe streaming approach, so pasting
// doesn't require holding every sample's decompressed column in memory
// at once (matters at panel scale, many samples).
//
// Same container/tool set as BINCOV_SET_BINS/BINCOV_MAKE_COLUMNS -- see
// those modules' notes.
//
process BINCOV_ZPASTE {
    tag "${meta.id}"
    label 'process_single'

    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/0b/0b4d52ca9a56d07be3f78a12af654e5116f5112908dba277e6796fd9dfb83fe5/data' :
        'community.wave.seqera.io/library/bcftools_htslib:1.23.1--9f08ec665533d64a' }"

    input:
    tuple val(meta), path(bin_locs), path(bincov_columns)

    output:
    tuple val(meta), path("*.RD.txt.gz"), path("*.RD.txt.gz.tbi"), emit: bincov_matrix

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p column_file_fifos
    FILE_NUM=0
    for COLUMN_FILE in ${bin_locs} ${bincov_columns}; do
        FIFO=\$(printf "column_file_fifos/%08d" \$FILE_NUM)
        mkfifo "\$FIFO"
        bgzip -@\$(nproc) -cd "\$COLUMN_FILE" > "\$FIFO" &
        FILE_NUM=\$((FILE_NUM + 1))
    done

    paste column_file_fifos/* | bgzip -@\$(nproc) -c > ${prefix}.RD.txt.gz
    tabix -p bed ${prefix}.RD.txt.gz
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo "" | gzip > ${prefix}.RD.txt.gz
    touch ${prefix}.RD.txt.gz.tbi
    """
}
