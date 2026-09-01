#!/usr/bin/env nextflow

//
// Standalone entry point: sample sheet -> per-sample PE/SR/SD/RD evidence
// files (see subworkflows/local/bam_collect_evidence).
//
// Run directly, e.g.:
//   nextflow run workflows/collect_evidence.nf -profile apptainer \
//     --input samplesheet.tsv --fasta ref.fa --fasta_fai ref.fa.fai \
//     --reference_dict ref.dict \
//     --sd_locs_vcf dbsnp.vcf --sd_locs_vcf_idx dbsnp.vcf.idx \
//     --preprocessed_intervals preprocessed_intervals.interval_list \
//     --outdir results
//
// This is the evidence-collection sibling to call_manta.nf/call_wham.nf --
// same tier (sample sheet -> per-sample output), not downstream of
// clustering/harmonisation. See
// ../README.md#composability and bam_collect_evidence/main.nf's own
// top-of-file note for why this exists ahead of genotyping.
//

include { UTILS_INPUT_CHANNELS  } from '../subworkflows/local/utils_input_channels/main'
include { BAM_COLLECT_EVIDENCE  } from '../subworkflows/local/bam_collect_evidence/main'

workflow {
    main:
    def pipelineRoot = "${projectDir}/.."

    [
        'reference_dict', 'sd_locs_vcf', 'sd_locs_vcf_idx', 'preprocessed_intervals'
    ].each { p ->
        if (!params[p]) {
            error "collect_evidence.nf requires --${p}"
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

    BAM_COLLECT_EVIDENCE.out.pe.view { meta, pe, tbi -> "PE evidence for ${meta.id}: ${pe}" }
    BAM_COLLECT_EVIDENCE.out.sr.view { meta, sr, tbi -> "SR evidence for ${meta.id}: ${sr}" }
    BAM_COLLECT_EVIDENCE.out.sd.view { meta, sd, tbi -> "SD evidence for ${meta.id}: ${sd}" }
    BAM_COLLECT_EVIDENCE.out.counts.view { meta, counts -> "RD counts for ${meta.id}: ${counts}" }
}
