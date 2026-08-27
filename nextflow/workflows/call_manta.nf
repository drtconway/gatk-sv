#!/usr/bin/env nextflow

//
// Standalone entry point: sample sheet -> Manta diploid SV VCFs.
//
// Run directly, e.g.:
//   nextflow run workflows/call_manta.nf -profile apptainer \
//     --input samplesheet.tsv --fasta ref.fa --fasta_fai ref.fa.fai
//
// This is the first slice of the composable pipeline described in
// ../README.md#composability: no panel bundle, harmonisation, or
// genotyping yet -- just the caller subworkflow wired to a sample sheet.
//

include { UTILS_INPUT_CHANNELS } from '../subworkflows/local/utils_input_channels/main'
include { BAM_CALL_MANTA       } from '../subworkflows/local/bam_call_manta/main'

workflow {
    main:
    // projectDir resolves to the directory of *this* script (workflows/)
    // when run standalone, so the pipeline root -- where assets/ lives --
    // has to be derived rather than assumed to be projectDir itself.
    def pipelineRoot = "${projectDir}/.."

    UTILS_INPUT_CHANNELS(params.input, pipelineRoot, params.fasta, params.fasta_fai)

    BAM_CALL_MANTA(
        UTILS_INPUT_CHANNELS.out.samples,
        UTILS_INPUT_CHANNELS.out.fasta,
        UTILS_INPUT_CHANNELS.out.fasta_fai
    )

    BAM_CALL_MANTA.out.vcf.view { meta, vcf -> "Manta diploid SV VCF for ${meta.id}: ${vcf}" }
}
