#!/usr/bin/env nextflow

//
// Standalone entry point: sample sheet -> standardized Manta SV VCFs.
//
// Run directly, e.g.:
//   nextflow run workflows/call_manta.nf -profile apptainer \
//     --input samplesheet.tsv --fasta ref.fa --fasta_fai ref.fa.fai \
//     --primary_contigs_fai contigs.fai \
//     --manta_region_bed manta_region.bed.gz \
//     --manta_region_bed_tbi manta_region.bed.gz.tbi
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

    ['primary_contigs_fai', 'manta_region_bed', 'manta_region_bed_tbi'].each { p ->
        if (!params[p]) {
            error "call_manta.nf requires --${p}"
        }
    }

    UTILS_INPUT_CHANNELS(params.input, pipelineRoot, params.fasta, params.fasta_fai)

    contigs_fai = channel.of([ [id: 'contigs'], file(params.primary_contigs_fai) ])
    manta_region_bed = channel.of([
        [id: 'manta_region'],
        file(params.manta_region_bed),
        file(params.manta_region_bed_tbi)
    ])

    BAM_CALL_MANTA(
        UTILS_INPUT_CHANNELS.out.samples,
        UTILS_INPUT_CHANNELS.out.fasta,
        UTILS_INPUT_CHANNELS.out.fasta_fai,
        manta_region_bed,
        contigs_fai,
        params.min_svsize
    )

    BAM_CALL_MANTA.out.vcf.view { meta, vcf -> "Standardized Manta VCF for ${meta.id}: ${vcf}" }
}
