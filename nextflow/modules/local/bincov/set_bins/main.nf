//
// Determine the common bin size/genomic bin grid from one sample's
// CollectReadCounts output (any one sample is representative -- all
// samples in a run were collected against the same
// params.preprocessed_intervals), for the panel-wide RD/bincov matrix's
// column-alignment reference. Reproduces wdl/MakeBincovMatrix.wdl's
// SetBins task, minus its bincov_matrix_samples/bincov_matrix handling
// (that exists to merge a *new* sample's counts against an *existing*
// panel matrix -- not something this pipeline does yet; revisit once a
// genotype-new-sample-style incremental path exists).
//
// Plain shell (awk/sed/bgzip), not a GATK walker -- same container as
// modules/local/tabix/modules/local/bgzip (bcftools/htslib), verified
// directly to have every tool this script uses.
//
process BINCOV_SET_BINS {
    tag "${meta.id}"
    label 'process_single'

    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/0b/0b4d52ca9a56d07be3f78a12af654e5116f5112908dba277e6796fd9dfb83fe5/data' :
        'community.wave.seqera.io/library/bcftools_htslib:1.23.1--9f08ec665533d64a' }"

    input:
    tuple val(meta), path(count_file)

    output:
    tuple val(meta), path("locs.bed.gz"), emit: bin_locs
    env 'BINSIZE'                       , emit: binsize

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    firstchar=\$(zcat -f ${count_file} | head -c 1)
    if [ "\$firstchar" == '@' ]; then
        shift=1
    else
        shift=0
    fi

    zcat -f ${count_file} \\
        | sed '/^@/d' \\
        | sed '/^CONTIG\tSTART\tEND\tCOUNT\$/d' \\
        | sed '/^#/d' \\
        | awk -v x="\${shift}" 'BEGIN{OFS="\\t"}{\$2=\$2-x; print \$1,\$2,\$3}' > tmp_locs

    # Use the most common interval size from the first 1000 bins, same
    # heuristic as wdl/MakeBincovMatrix.wdl's SetBins (no fixed binsize
    # override wired through here -- GATK-SV's own binsize param is
    # optional and this pipeline doesn't currently need to pin one).
    export BINSIZE=\$(sed -n '1,1000p' tmp_locs | awk '{ print \$3-\$2 }' | sort | uniq -c | sort -nrk1,1 | sed -n '1p' | awk '{ print \$2 }')

    awk -v FS='\\t' -v b="\${BINSIZE}" 'BEGIN{ print "#Chr\\tStart\\tEnd" } { if (\$3-\$2==b) print \$0 }' tmp_locs \\
        | bgzip -c > locs.bed.gz
    """

    stub:
    """
    export BINSIZE=1000
    echo "" | gzip > locs.bed.gz
    """
}
