//
// Give every interval the full sample list (to genotype everyone), then
// split into small (<5kb) and large intervals, in chunks of n_per_split
// lines each. Reproduces wdl/GenotypeCpxCnvsPerBatch.wdl's SplitBedBySize,
// with the sample list taken from the median coverage file's header
// (GetSampleIdsFromMedianCoverageFile) in the same task.
//
// Always emits at least lt5kb.1 and gt5kb.1, possibly empty, as upstream;
// RdTest handles an empty chunk.
//
process SPLIT_BED_BY_SIZE {
    tag "${meta.id}"
    label 'process_single'

    container 'us.gcr.io/broad-dsde-methods/gatk-sv/sv-pipeline:2026-06-02-v1.1-24f0fd08'

    input:
    tuple val(meta), path(bed)
    path(median_file)
    val(n_per_split_small)
    val(n_per_split_large)

    output:
    tuple val(meta), path("{lt,gt}5kb.*"), emit: beds

    when:
    task.ext.when == null || task.ext.when

    script:
    // Size cutoff exactly as upstream: size_cutoff_kb (5) is compared with
    // the length in bp, never multiplied by 1000, so in practice every
    // interval goes to gt5kb (chunks of n_per_split_large). Kept as is to
    // match upstream's output.
    """
    head -1 ${median_file} | sed -e 's/\\t/\\n/g' > samples.list

    zcat ${bed} \\
        | awk -F "\\t" -v OFS="\\t" -v samples=\$(paste -s -d, samples.list) \\
            '{ print \$1, \$2, \$3, \$4, samples, \$6 }' \\
        | sort -T . -Vk1,1 -k2,2n -k3,3n \\
        > merged.bed

    touch lt5kb.bed gt5kb.bed
    awk -F "\\t" -v OFS="\\t" '{ if (\$3 - \$2 < 5) { print > "lt5kb.bed" } else { print > "gt5kb.bed" } }' merged.bed
    rm -f merged.bed

    printf 'lt5kb.bed %s\\ngt5kb.bed %s\\n' ${n_per_split_small} ${n_per_split_large} \\
        | while read SIZED_BED LINES_PER_SPLIT; do
            PREFIX="\${SIZED_BED%.*}."
            N_LINES=\$(wc -l < \$SIZED_BED)
            N_CHUNKS=\$((N_LINES / LINES_PER_SPLIT))
            if [ \$N_CHUNKS -gt 1 ]; then
                N_DIGITS=\${#N_CHUNKS}
                split -d -a \$N_DIGITS -n l/\$N_CHUNKS \\
                    --numeric-suffixes=\$(printf "%0\${N_DIGITS}d" 1) \\
                    \$SIZED_BED \\
                    \$PREFIX
                rm -f \$SIZED_BED
            else
                mv \$SIZED_BED \${PREFIX}1
            fi
        done
    """

    stub:
    """
    touch lt5kb.1 gt5kb.1
    """
}
