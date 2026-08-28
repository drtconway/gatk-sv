#!/usr/bin/env nextflow

//
// Standalone entry point: sample sheet -> standardized Wham SV VCFs.
//
// Run directly, e.g.:
//   nextflow run workflows/call_wham.nf -profile apptainer \
//     --input samplesheet.tsv --fasta ref.fa --fasta_fai ref.fa.fai \
//     --primary_contigs_list primary_contigs.list \
//     --primary_contigs_fai contigs.fai
//
// Same composable-caller pattern as call_manta.nf -- see
// ../README.md#composability.
//

include { UTILS_INPUT_CHANNELS } from '../subworkflows/local/utils_input_channels/main'
include { BAM_CALL_WHAM        } from '../subworkflows/local/bam_call_wham/main'

workflow {
    main:
    // projectDir resolves to the directory of *this* script (workflows/)
    // when run standalone, so the pipeline root -- where assets/ lives --
    // has to be derived rather than assumed to be projectDir itself.
    def pipelineRoot = "${projectDir}/.."

    // primary_contigs_list (plain list, for whamg -c) and
    // primary_contigs_fai (.fai format, for svtk standardize --contigs)
    // are both needed and serve different tools -- see nextflow.config.
    ['primary_contigs_list', 'primary_contigs_fai'].each { p ->
        if (!params[p]) {
            error "call_wham.nf requires --${p}"
        }
    }

    UTILS_INPUT_CHANNELS(params.input, pipelineRoot, params.fasta, params.fasta_fai)

    contigs_fai = channel.of([ [id: 'contigs'], file(params.primary_contigs_fai) ])

    BAM_CALL_WHAM(
        UTILS_INPUT_CHANNELS.out.samples,
        UTILS_INPUT_CHANNELS.out.fasta,
        UTILS_INPUT_CHANNELS.out.fasta_fai,
        contigs_fai,
        params.min_svsize
    )

    BAM_CALL_WHAM.out.vcf.view { meta, vcf -> "Standardized Wham VCF for ${meta.id}: ${vcf}" }
}
