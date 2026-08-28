//
// Shared input parsing for standalone workflows/call_*.nf entry points:
// sample sheet -> per-sample [meta, bam, bai] channel, plus the reference
// genome as broadcastable [meta, path] channels.
//
// Factored out once a second caller (Wham) needed the exact same ~10
// lines call_manta.nf already had -- see README.md#composability.
//

include { samplesheetToList } from 'plugin/nf-schema'

workflow UTILS_INPUT_CHANNELS {
    take:
    samplesheet_path  // val: params.input
    pipeline_root      // val: absolute path to the nextflow/ pipeline root (see note below)
    fasta_path          // val: params.fasta
    fasta_fai_path      // val: params.fasta_fai

    main:
    // 'family_id', 'sample_id', and 'ped_id' are all tagged meta: in the
    // schema, so samplesheetToList already merges them into one meta map
    // per row: [ [family_id:.., id:.., ped_id:..], bam, bai ]
    //
    // ped_id is optional and defaults to sample_id (meta.id) when absent.
    // The two commonly differ: sample_id is often a lab/external
    // reference (e.g. "23PS00385-RDN0206"), while a pedigree file may use
    // a different internal family/relationship-position convention (e.g.
    // "RDN0206-00"). Every PED-keyed lookup downstream (ploidy tables,
    // VCF sample-column rewrites via bcftools reheader) must use ped_id,
    // not sample_id -- a real KeyError from exactly this mismatch is what
    // surfaced the need for this field. sample_id remains meta.id and is
    // used for filenames/display; ped_id is threaded through separately
    // wherever PED-derived data is consumed.
    samples = channel
        .fromList(samplesheetToList(samplesheet_path, "${pipeline_root}/assets/schema_samplesheet.json"))
        .map { meta, bam, bai ->
            def meta_with_ped_id = meta.ped_id ? meta : meta + [ ped_id: meta.id ]
            [ meta_with_ped_id, file(bam), file(bai) ]
        }

    fasta     = channel.of([ [ id: 'reference' ], file(fasta_path) ])
    fasta_fai = channel.of([ [ id: 'reference' ], file(fasta_fai_path) ])

    emit:
    samples    // channel: [ meta, bam, bai ]
    fasta      // channel: [ meta, fasta ]
    fasta_fai  // channel: [ meta, fasta_fai ]
}
