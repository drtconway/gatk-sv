//
// Manta germline SV calling, one BAM/CRAM per sample.
//
// GATK-SV's own Manta task (wdl/Manta.wdl) runs one Manta invocation per
// sample and keeps only the diploid SV VCF -- candidate VCFs are Manta
// intermediates GATK-SV never consumes downstream, so we drop them here
// too rather than plumbing them through for no consumer.
//

include { MANTA_GERMLINE } from '../../../modules/nf-core/manta/germline/main'

workflow BAM_CALL_MANTA {
    take:
    bam        // channel: [mandatory] [ meta, bam, bai ]
    fasta      // channel: [mandatory] [ meta2, fasta ]
    fasta_fai  // channel: [mandatory] [ meta3, fai ]

    main:
    // Manta's target_bed/target_bed_tbi are for restricting calling to
    // intervals (e.g. exome). GATK-SV runs Manta genome-wide, so we pass
    // no target bed.
    manta_input = bam.map { meta, bam_file, bai_file ->
        [ meta, bam_file, bai_file, [], [] ]
    }

    // fasta/fasta_fai are one reference shared across every sample, but
    // MANTA_GERMLINE declares them as ordinary (non-value) channel inputs.
    // Without .collect() they're each consumed once and then exhausted, so
    // Nextflow silently pairs positionally and drops every sample after
    // the first -- .collect() turns them into a single broadcastable value
    // reused for every item in manta_input.
    MANTA_GERMLINE(
        manta_input,
        fasta.collect(),
        fasta_fai.collect(),
        []  // no config override
    )
    // NB: MANTA_GERMLINE.out.versions_manta is deliberately not consumed.
    // It's a `topic: versions, emit: ..., eval("configManta.py --version")`
    // output -- Nextflow's `eval` output type always shells out for real,
    // even under -stub-run, so referencing it breaks stub-mode validation
    // when the real tool isn't installed. Wire it up once we're running
    // against real containers and actually need version tracking.

    vcf = MANTA_GERMLINE.out.diploid_sv_vcf
    vcf_tbi = MANTA_GERMLINE.out.diploid_sv_vcf_tbi

    emit:
    vcf       // channel: [ meta, vcf ]
    vcf_tbi   // channel: [ meta, vcf_tbi ]
}
