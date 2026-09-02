//
// Merge every sample's RD (binned read-depth) counts, from
// bam_collect_evidence, into one panel-wide bincov matrix (one column
// per sample, sharing a common bin grid). Mirrors GATK-SV's
// MakeBincovMatrix workflow (wdl/MakeBincovMatrix.wdl) -- see
// vcfs_merge_evidence/main.nf's own note for why this is a *separate*
// subworkflow rather than folded into that one: this is plain shell
// (awk/paste/bgzip), not a GATK walker, and needs a "pick one sample as
// the reference grid, then fan out" shape none of this pipeline's other
// merge steps have.
//
// Deliberately not implemented: the WDL's bincov_matrix/
// bincov_matrix_samples inputs (merging a *new* sample's counts against
// an already-existing panel matrix, for incremental runs) -- this
// subworkflow only builds a matrix from scratch across whatever samples
// are in the current run. Revisit once a genotype-new-sample-style
// incremental path exists (see README.md's panel-based design).
//

include { BINCOV_SET_BINS     } from '../../../modules/local/bincov/set_bins/main'
include { BINCOV_MAKE_COLUMNS } from '../../../modules/local/bincov/make_columns/main'
include { BINCOV_ZPASTE       } from '../../../modules/local/bincov/zpaste/main'

workflow VCFS_MERGE_READ_COUNTS {
    take:
    counts        // channel: [ meta, counts_tsv_gz ] -- one per sample, from bam_collect_evidence
    cohort_name    // val: String, used as the output prefix

    main:
    // Any one sample's counts are representative for establishing the
    // bin grid -- every sample in a run was collected against the same
    // params.preprocessed_intervals (see bam_collect_evidence). `.first()`
    // takes the first element of the collected channel deterministically
    // (Nextflow channels preserve arrival order for a single upstream
    // process here, since GATK_COLLECT_READ_COUNTS has no scatter of its
    // own per sample beyond the one call).
    reference_counts = counts.first()

    BINCOV_SET_BINS(reference_counts)
    // Both broadcast across every BINCOV_MAKE_COLUMNS call below, same
    // .collect()-as-broadcast pattern used throughout this pipeline for
    // a single shared value (see e.g. bam_call_wham's fasta/fasta_fai) --
    // BINCOV_SET_BINS runs exactly once (on reference_counts), so its
    // output channels each carry exactly one value to broadcast.
    bin_locs = BINCOV_SET_BINS.out.bin_locs.map { meta, locs -> locs }.collect()
    // .collect() on a single-value channel gives a one-element list --
    // unwrap it back to a plain value for the `val(binsize)` process
    // input below (path inputs tolerate a collected list directly, but a
    // val input needs the bare value, not [value]).
    binsize  = BINCOV_SET_BINS.out.binsize.collect().map { it[0] }

    BINCOV_MAKE_COLUMNS(
        counts,
        binsize,
        bin_locs
    )

    // .combine() flattens list-valued items together with whatever it's
    // combined against, rather than nesting them -- e.g.
    // channel.of([1,2]).combine(channel.of('x')) emits [1,2,'x'], not
    // [[1,2],'x']. Wrapping the columns list in an extra list first
    // (.map { l -> [l] }) is what keeps it as one distinguishable tuple
    // element after combine, confirmed empirically (not documented
    // clearly) before relying on it here.
    bin_locs_value = bin_locs.map { it[0] }

    zpaste_input = BINCOV_MAKE_COLUMNS.out.bincov_column
        .map { meta, col -> col }
        .collect()
        .map { cols -> [ cols ] }
        .combine(bin_locs_value)
        .map { cols, locs -> [ [id: cohort_name], locs, cols ] }

    BINCOV_ZPASTE(zpaste_input)

    merged_counts = BINCOV_ZPASTE.out.bincov_matrix

    emit:
    merged_counts  // channel: [ meta, RD_txt_gz, RD_txt_gz_tbi ]
}
