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
    // 'family_id' and 'sample_id' are both tagged meta: in the schema, so
    // samplesheetToList already merges them into one meta map per row:
    // [ [family_id:.., id:..], bam, bai ]
    samples = channel
        .fromList(samplesheetToList(samplesheet_path, "${pipeline_root}/assets/schema_samplesheet.json"))
        .map { meta, bam, bai -> [ meta, file(bam), file(bai) ] }

    fasta     = channel.of([ [ id: 'reference' ], file(fasta_path) ])
    fasta_fai = channel.of([ [ id: 'reference' ], file(fasta_fai_path) ])

    emit:
    samples    // channel: [ meta, bam, bai ]
    fasta      // channel: [ meta, fasta ]
    fasta_fai  // channel: [ meta, fasta_fai ]
}
