//
// Merge every sample's PE/SR/SD evidence (from bam_collect_evidence)
// across a run into one panel-wide PE file, one panel-wide SR file, and
// one panel-wide BAF file. Mirrors GATK-SV's BatchEvidenceMerging
// workflow (wdl/BatchEvidenceMerging.wdl) -- see
// README.md#status for where this fits (the next unbuilt prerequisite
// for genotyping, after per-sample evidence collection).
//
// Deliberately simpler than the WDL in one respect: GATK-SV's
// MergeEvidence/SDtoBAF tasks both have an optional rename_samples step
// (an awk rewrite of each file's sample column to the batch's canonical
// sample list) because their per-sample evidence files are keyed by
// whatever sample_id was passed to GatherSampleEvidence, which may not
// match the batch's own naming. We don't have that mismatch to correct:
// bam_collect_evidence's GATK_COLLECT_SV_EVIDENCE already writes every
// PE/SR/SD file with --sample-name set to meta.ped_id at collection time
// (the same PED-keyed identity every other lookup in this pipeline
// uses), so there's nothing left to rename here. Also not implemented:
// GATK-SV's subset_primary_contigs option (drops ALT-contig records
// before merging) -- not exercised yet since evidence collection here
// isn't restricted to primary contigs either; revisit together if it
// becomes necessary, same as the primary-contigs handling noted
// elsewhere in this pipeline (see README.md's ALT-contig bug entries).
//
// BAF, not raw BAF files: GATK-SV's BatchEvidenceMerging accepts either
// pre-existing BAF_files or SD_files (converted via SDtoBAF) -- we only
// ever produce SD (bam_collect_evidence has no BAF collection step), so
// SDtoBAF (GATK_SITE_DEPTH_TO_BAF) is the only path implemented, not the
// PrintSVEvidence-with-evidence="baf" alternative the WDL also supports.
//

include { GATK_PRINT_SV_EVIDENCE as GATK_PRINT_SV_EVIDENCE_PE } from '../../../modules/local/gatk4/print_sv_evidence/main'
include { GATK_PRINT_SV_EVIDENCE as GATK_PRINT_SV_EVIDENCE_SR } from '../../../modules/local/gatk4/print_sv_evidence/main'
include { GATK_SITE_DEPTH_TO_BAF } from '../../../modules/local/gatk4/site_depth_to_baf/main'
include { TABIX_SV_EVIDENCE as TABIX_MERGED_PE  } from '../../../modules/local/tabix_sv_evidence/main'
include { TABIX_SV_EVIDENCE as TABIX_MERGED_SR  } from '../../../modules/local/tabix_sv_evidence/main'
include { TABIX_SV_EVIDENCE as TABIX_MERGED_BAF } from '../../../modules/local/tabix_sv_evidence/main'

workflow VCFS_MERGE_EVIDENCE {
    take:
    pe               // channel: [ meta, pe_txt_gz, pe_txt_gz_tbi ] -- one per sample, from bam_collect_evidence
    sr               // channel: [ meta, sr_txt_gz, sr_txt_gz_tbi ] -- one per sample
    sd               // channel: [ meta, sd_txt_gz, sd_txt_gz_tbi ] -- one per sample
    dict             // channel: [ meta2, reference_dict ]
    sd_locs_vcf      // channel: [ meta3, vcf, vcf_idx ] -- e.g. params.sd_locs_vcf (dbSNP sites);
                      // NB vcf_idx is a GATK .idx sibling, not a bgzip .tbi -- see
                      // bam_collect_evidence/main.nf's own note.
    cohort_name       // val: String, used as the output prefix

    main:
    // Each merge collects across all samples into one call -- same
    // "single genome-wide call, not scattered by contig" simplification
    // already made throughout this pipeline (see e.g.
    // vcfs_cluster_svcluster's own note on SVCluster).
    pe_input = pe
        .map { meta, f, tbi -> [ f, tbi ] }
        .collect(flat: false)
        .map { pairs -> [ [id: "${cohort_name}.pe"], pairs.collect { it[0] }, pairs.collect { it[1] } ] }

    GATK_PRINT_SV_EVIDENCE_PE(pe_input, dict.map { meta, d -> d }.collect())
    TABIX_MERGED_PE(GATK_PRINT_SV_EVIDENCE_PE.out.evidence)

    sr_input = sr
        .map { meta, f, tbi -> [ f, tbi ] }
        .collect(flat: false)
        .map { pairs -> [ [id: "${cohort_name}.sr"], pairs.collect { it[0] }, pairs.collect { it[1] } ] }

    GATK_PRINT_SV_EVIDENCE_SR(sr_input, dict.map { meta, d -> d }.collect())
    TABIX_MERGED_SR(GATK_PRINT_SV_EVIDENCE_SR.out.evidence)

    // Plain cohort_name (not "${cohort_name}.baf" like the pe/sr meta.id
    // choices above): GATK_SITE_DEPTH_TO_BAF's own script already appends
    // ".baf.txt.gz" to its prefix (mirroring wdl/BatchEvidenceMerging.wdl's
    // SDtoBAF: "~{batch}.baf.txt.gz") -- adding ".baf" here too produced a
    // real "cohort.baf.baf.txt.gz" double-suffix bug, caught by
    // -stub-run's own file listing, not a wiring failure (every process
    // still succeeded; the output filename was just wrong).
    sd_input = sd
        .map { meta, f, tbi -> [ f, tbi ] }
        .collect(flat: false)
        .map { pairs -> [ [id: cohort_name], pairs.collect { it[0] }, pairs.collect { it[1] } ] }

    GATK_SITE_DEPTH_TO_BAF(
        sd_input,
        dict.map { meta, d -> d }.collect(),
        sd_locs_vcf.map { meta, v, idx -> v }.collect(),
        sd_locs_vcf.map { meta, v, idx -> idx }.collect()
    )
    TABIX_MERGED_BAF(GATK_SITE_DEPTH_TO_BAF.out.baf)

    merged_pe  = TABIX_MERGED_PE.out.evidence
    merged_sr  = TABIX_MERGED_SR.out.evidence
    merged_baf = TABIX_MERGED_BAF.out.evidence

    emit:
    merged_pe   // channel: [ meta, pe_txt_gz, pe_txt_gz_tbi ]
    merged_sr   // channel: [ meta, sr_txt_gz, sr_txt_gz_tbi ]
    merged_baf  // channel: [ meta, baf_txt_gz, baf_txt_gz_tbi ]
}
