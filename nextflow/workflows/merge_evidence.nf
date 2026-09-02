#!/usr/bin/env nextflow

//
// Standalone entry point: sample sheet -> panel-wide merged PE/SR/BAF/RD
// evidence (see subworkflows/local/vcfs_merge_evidence and
// subworkflows/local/vcfs_merge_read_counts).
//
// Run directly, e.g.:
//   nextflow run workflows/merge_evidence.nf -profile apptainer \
//     --input samplesheet.tsv --fasta ref.fa --fasta_fai ref.fa.fai \
//     --reference_dict ref.dict \
//     --sd_locs_vcf dbsnp.vcf --sd_locs_vcf_idx dbsnp.vcf.idx \
//     --preprocessed_intervals preprocessed_intervals.interval_list \
//     --outdir results
//
// Composes bam_collect_evidence (per-sample) with vcfs_merge_evidence +
// vcfs_merge_read_counts (panel-wide) -- same "one entry point runs the
// full chain" pattern as combine_batches.nf. This is the last unbuilt
// prerequisite for genotyping (TrainSVGenotyping/GenotypeSVs, not yet
// implemented) -- see README.md#status for the full reasoning chain.
//

include { UTILS_INPUT_CHANNELS     } from '../subworkflows/local/utils_input_channels/main'
include { BAM_COLLECT_EVIDENCE     } from '../subworkflows/local/bam_collect_evidence/main'
include { VCFS_MERGE_EVIDENCE      } from '../subworkflows/local/vcfs_merge_evidence/main'
include { VCFS_MERGE_READ_COUNTS   } from '../subworkflows/local/vcfs_merge_read_counts/main'

workflow {
    main:
    def pipelineRoot = "${projectDir}/.."

    [
        'reference_dict', 'sd_locs_vcf', 'sd_locs_vcf_idx', 'preprocessed_intervals'
    ].each { p ->
        if (!params[p]) {
            error "merge_evidence.nf requires --${p}"
        }
    }

    UTILS_INPUT_CHANNELS(params.input, pipelineRoot, params.fasta, params.fasta_fai)

    dict = channel.of([ [id: 'reference'], file(params.reference_dict) ])
    // NB sd_locs_vcf_idx is a GATK .idx sibling, not a bgzip .tbi -- see
    // bam_collect_evidence/main.nf's own note on this param.
    sd_locs_vcf = channel.of([
        [id: 'sd_locs'],
        file(params.sd_locs_vcf),
        file(params.sd_locs_vcf_idx)
    ])
    preprocessed_intervals = channel.of([ [id: 'intervals'], file(params.preprocessed_intervals) ])

    BAM_COLLECT_EVIDENCE(
        UTILS_INPUT_CHANNELS.out.samples,
        UTILS_INPUT_CHANNELS.out.fasta,
        UTILS_INPUT_CHANNELS.out.fasta_fai,
        dict,
        sd_locs_vcf,
        preprocessed_intervals
    )

    VCFS_MERGE_EVIDENCE(
        BAM_COLLECT_EVIDENCE.out.pe,
        BAM_COLLECT_EVIDENCE.out.sr,
        BAM_COLLECT_EVIDENCE.out.sd,
        dict,
        sd_locs_vcf,
        'cohort'
    )

    VCFS_MERGE_READ_COUNTS(
        BAM_COLLECT_EVIDENCE.out.counts,
        'cohort'
    )

    VCFS_MERGE_EVIDENCE.out.merged_pe.view { meta, pe, tbi -> "Merged PE evidence: ${pe}" }
    VCFS_MERGE_EVIDENCE.out.merged_sr.view { meta, sr, tbi -> "Merged SR evidence: ${sr}" }
    VCFS_MERGE_EVIDENCE.out.merged_baf.view { meta, baf, tbi -> "Merged BAF evidence: ${baf}" }
    VCFS_MERGE_READ_COUNTS.out.merged_counts.view { meta, rd, tbi -> "Merged RD bincov matrix: ${rd}" }
}
