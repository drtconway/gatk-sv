//
// Per-sample PE/SR/SD/RD evidence collection from a BAM/CRAM. This is the
// evidence-collection sibling to bam_call_manta/bam_call_wham -- same
// tier (one BAM in, per-sample output out), run alongside them rather
// than downstream of clustering/harmonisation. Mirrors GATK-SV's
// GatherSampleEvidence workflow (wdl/GatherSampleEvidence.wdl), narrowed
// to just the SV-evidence + coverage collection it does (not the SV
// calling it also bundles -- that's already covered by
// bam_call_manta/bam_call_wham here).
//
// This evidence (raw per-sample PE/SR/SD/RD) is what a future genotyping
// stage (GATK's TrainSVGenotyping/GenotypeSVs walkers, not yet
// implemented) needs, once merged panel-wide -- see
// README.md#not-yet-started. Collecting it per-sample now, ahead of
// building the merge step or genotyping itself, follows the same
// incremental "build one piece, validate on real HPC data, then build
// the next" rhythm the rest of this pipeline has used.
//
// Two of GATK's walkers, chained with small format-conversion steps
// because neither the GATK image nor a single "do everything" container
// covers what both need:
//   1. CollectSVEvidence (GATK_COLLECT_SV_EVIDENCE) -> PE/SR/SD, already
//      bgzipped by GATK's own writer, indexed here (TABIX_SV_EVIDENCE)
//      since `tabix` isn't in the GATK image.
//   2. CollectReadCounts (GATK_COLLECT_READ_COUNTS) -> RD (binned
//      coverage), plain TSV, bgzipped here (BGZIP) since `bgzip` isn't in
//      the GATK image either. Not indexed -- GATK-SV's own downstream
//      evidence-merging step consumes this as plain bgzipped TSV, not a
//      tabix-indexed file (see wdl/GatherBatchEvidence.wdl's
//      BatchEvidenceMerging).
//

include { GATK_COLLECT_SV_EVIDENCE } from '../../../modules/local/gatk4/collect_sv_evidence/main'
include { GATK_COLLECT_READ_COUNTS } from '../../../modules/local/gatk4/collect_read_counts/main'
include { TABIX_SV_EVIDENCE as TABIX_PE } from '../../../modules/local/tabix_sv_evidence/main'
include { TABIX_SV_EVIDENCE as TABIX_SR } from '../../../modules/local/tabix_sv_evidence/main'
include { TABIX_SV_EVIDENCE as TABIX_SD } from '../../../modules/local/tabix_sv_evidence/main'
include { BGZIP as BGZIP_COUNTS } from '../../../modules/local/bgzip/main'

workflow BAM_COLLECT_EVIDENCE {
    take:
    bam                    // channel: [mandatory] [ meta, bam, bai ]
    fasta                  // channel: [mandatory] [ meta2, fasta ]
    fasta_fai              // channel: [mandatory] [ meta3, fasta_fai ]
    dict                    // channel: [mandatory] [ meta4, reference_dict ]
    sd_locs_vcf              // channel: [mandatory] [ meta5, vcf, vcf_idx ] -- e.g. params.sd_locs_vcf (dbSNP sites);
                              // NB vcf_idx is a GATK .idx sibling, not a bgzip .tbi -- see
                              // modules/local/gatk4/collect_sv_evidence's own note.
    preprocessed_intervals    // channel: [mandatory] [ meta6, interval_list ]

    main:
    // fasta/fasta_fai/dict/sd_locs_vcf/preprocessed_intervals are all one
    // value shared across every sample -- .collect() broadcasts them the
    // same way every other per-sample subworkflow in this pipeline does
    // (see bam_call_manta/bam_call_wham for the same pattern).
    GATK_COLLECT_SV_EVIDENCE(
        bam,
        fasta.map { meta, f -> f }.collect(),
        fasta_fai.map { meta, f -> f }.collect(),
        dict.map { meta, d -> d }.collect(),
        sd_locs_vcf.map { meta, v, idx -> v }.collect(),
        sd_locs_vcf.map { meta, v, idx -> idx }.collect()
    )

    TABIX_PE(GATK_COLLECT_SV_EVIDENCE.out.pe)
    TABIX_SR(GATK_COLLECT_SV_EVIDENCE.out.sr)
    TABIX_SD(GATK_COLLECT_SV_EVIDENCE.out.sd)

    GATK_COLLECT_READ_COUNTS(
        bam,
        fasta.map { meta, f -> f }.collect(),
        fasta_fai.map { meta, f -> f }.collect(),
        dict.map { meta, d -> d }.collect(),
        preprocessed_intervals.map { meta, i -> i }.collect()
    )

    BGZIP_COUNTS(GATK_COLLECT_READ_COUNTS.out.counts)

    pe    = TABIX_PE.out.evidence
    sr    = TABIX_SR.out.evidence
    sd    = TABIX_SD.out.evidence
    counts = BGZIP_COUNTS.out.gz

    emit:
    pe      // channel: [ meta, pe_txt_gz, pe_txt_gz_tbi ]
    sr      // channel: [ meta, sr_txt_gz, sr_txt_gz_tbi ]
    sd      // channel: [ meta, sd_txt_gz, sd_txt_gz_tbi ]
    counts  // channel: [ meta, counts_tsv_gz ]
}
