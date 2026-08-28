#!/usr/bin/env nextflow

//
// Standalone entry point: sample sheet -> Wham VCFs (contig-restricted,
// sample ID/TAGS fixed to match the pipeline's sample_id).
//
// Run directly, e.g.:
//   nextflow run workflows/call_wham.nf -profile apptainer \
//     --input samplesheet.tsv --fasta ref.fa --fasta_fai ref.fa.fai \
//     --primary_contigs_list primary_contigs.list
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

    if (!params.primary_contigs_list) {
        error "call_wham.nf requires --primary_contigs_list (see conf/modules.config)"
    }

    UTILS_INPUT_CHANNELS(params.input, pipelineRoot, params.fasta, params.fasta_fai)

    contig_list = channel.of([ [id: 'contigs'], file(params.primary_contigs_list) ])

    BAM_CALL_WHAM(
        UTILS_INPUT_CHANNELS.out.samples,
        UTILS_INPUT_CHANNELS.out.fasta,
        UTILS_INPUT_CHANNELS.out.fasta_fai,
        contig_list
    )

    BAM_CALL_WHAM.out.vcf.view { meta, vcf -> "Wham VCF for ${meta.id}: ${vcf}" }
}
