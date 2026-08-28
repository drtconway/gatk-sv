#!/usr/bin/env nextflow

//
// Standalone entry point: sample sheet -> VCFs from every implemented
// caller (currently Manta, Wham).
//
// Run directly, e.g.:
//   nextflow run workflows/call_all.nf -profile apptainer \
//     --input samplesheet.tsv --fasta ref.fa --fasta_fai ref.fa.fai \
//     --primary_contigs_list primary_contigs.list \
//     --manta_region_bed manta_region.bed.gz \
//     --manta_region_bed_tbi manta_region.bed.gz.tbi
//
// This is not build_panel.nf / genotype_new_sample.nf -- no panel bundle,
// harmonisation (SVCluster), or genotyping yet. It's the multi-caller
// layer underneath those: each caller subworkflow still runs
// independently and is still runnable standalone via its own
// workflows/call_<tool>.nf (see ../README.md#composability). This just
// avoids parsing the sample sheet and re-running every caller by hand
// when you want output from all of them for the same samples.
//
// Add a new `include` + call block here as each new caller subworkflow
// lands (Scramble, cn.MOPS, gCNV, ...).
//

include { UTILS_INPUT_CHANNELS } from '../subworkflows/local/utils_input_channels/main'
include { BAM_CALL_MANTA       } from '../subworkflows/local/bam_call_manta/main'
include { BAM_CALL_WHAM        } from '../subworkflows/local/bam_call_wham/main'

workflow {
    main:
    // projectDir resolves to the directory of *this* script (workflows/)
    // when run standalone, so the pipeline root -- where assets/ lives --
    // has to be derived rather than assumed to be projectDir itself.
    def pipelineRoot = "${projectDir}/.."

    ['primary_contigs_list', 'manta_region_bed', 'manta_region_bed_tbi'].each { p ->
        if (!params[p]) {
            error "call_all.nf requires --${p}"
        }
    }

    UTILS_INPUT_CHANNELS(params.input, pipelineRoot, params.fasta, params.fasta_fai)
    samples   = UTILS_INPUT_CHANNELS.out.samples
    fasta     = UTILS_INPUT_CHANNELS.out.fasta
    fasta_fai = UTILS_INPUT_CHANNELS.out.fasta_fai

    contig_list = channel.of([ [id: 'contigs'], file(params.primary_contigs_list) ])
    manta_region_bed = channel.of([
        [id: 'manta_region'],
        file(params.manta_region_bed),
        file(params.manta_region_bed_tbi)
    ])

    BAM_CALL_MANTA(samples, fasta, fasta_fai, manta_region_bed, contig_list)
    BAM_CALL_WHAM(samples, fasta, fasta_fai, contig_list)

    BAM_CALL_MANTA.out.vcf.view { meta, vcf -> "Manta diploid SV VCF for ${meta.id}: ${vcf}" }
    BAM_CALL_WHAM.out.vcf.view  { meta, vcf -> "Wham VCF for ${meta.id}: ${vcf}" }
}
