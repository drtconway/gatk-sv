//
// Genotype SVs at a sites-only VCF's sites, using panel-wide PE/SR/RD
// evidence + rf_cutoffs (from FilterBatchSites/AdjudicateSV) to train
// per-cohort genotyping cutoff tables, then apply them per contig.
// Mirrors GATK-SV's GenotypeBatch workflow (wdl/GenotypeBatch.wdl) --
// see README.md#path-to-genotyping-the-real-prerequisite-chain for the
// full chain this is the last link of.
//
// Chain: TrainSVGenotyping (once, whole VCF) -> per contig
// (PrintSVEvidence x3, subsetting panel-wide RD/PE/SR to that contig,
// then tabix-indexed -> GenotypeSVs, using the trained tables) -> concat
// contig shards back -> SeparateDepthPesr (split by INFO/ALGORITHMS) ->
// GenerateRegenoCoverageMedians (RD_MCR extraction, for RegenotypeCNVs --
// not yet built here, kept for parity with upstream's own output set).
//
// PrintSVEvidence's own per-contig output filenames must literally end
// in ".rd.txt.gz"/".pe.txt.gz"/".sr.txt.gz" (case-insensitive, but some
// prefix segment before that suffix is required) -- GATK's own codec
// detection is filename-suffix-based, confirmed empirically (a bare
// "rd.txt.gz" with no prefix fails "no suitable codecs found"; any real
// prefix works). Each of the three GATK_PRINT_SV_EVIDENCE_CONTIG aliases
// below gets its own task.ext.prefix in conf/modules.config accordingly.
//
// Every per-contig process below is keyed by the contig's own name
// (meta.id = contig), not the "${cohort_name}.genotype_batch.${contig}"
// scheme used elsewhere -- joining the three independent PrintSVEvidence
// scatters (RD/PE/SR) back together needs a shared key, and channel
// scatter/gather order across three separate process calls isn't
// guaranteed to line up positionally (same class of ordering risk
// vcfs_generate_batch_metrics's own CONCAT_OUTPUT_VCFS hit -- a `.join()`
// on a real key sidesteps it entirely rather than relying on order).
//

include { GATK_TRAINSVGENOTYPING } from '../../../modules/local/gatk4/train_sv_genotyping/main'
include { GATK_PRINT_SV_EVIDENCE_CONTIG as GATK_PRINT_SV_EVIDENCE_CONTIG_RD } from '../../../modules/local/gatk4/print_sv_evidence_contig/main'
include { GATK_PRINT_SV_EVIDENCE_CONTIG as GATK_PRINT_SV_EVIDENCE_CONTIG_PE } from '../../../modules/local/gatk4/print_sv_evidence_contig/main'
include { GATK_PRINT_SV_EVIDENCE_CONTIG as GATK_PRINT_SV_EVIDENCE_CONTIG_SR } from '../../../modules/local/gatk4/print_sv_evidence_contig/main'
include { TABIX_SV_EVIDENCE as TABIX_CONTIG_RD } from '../../../modules/local/tabix_sv_evidence/main'
include { TABIX_SV_EVIDENCE as TABIX_CONTIG_PE } from '../../../modules/local/tabix_sv_evidence/main'
include { TABIX_SV_EVIDENCE as TABIX_CONTIG_SR } from '../../../modules/local/tabix_sv_evidence/main'
include { GATK_GENOTYPESVS } from '../../../modules/local/gatk4/genotype_svs/main'
include { CONCAT_VCFS as CONCAT_GENOTYPE_BATCH_SHARDS } from '../../../modules/local/concat_vcfs/main'
include { SEPARATE_DEPTH_PESR } from '../../../modules/local/separate_depth_pesr/main'
include { GENERATE_REGENO_COVERAGE_MEDIANS } from '../../../modules/local/generate_regeno_coverage_medians/main'

workflow GENOTYPE_BATCH {
    take:
    vcf                 // channel: [ meta, vcf, vcf_tbi ] -- sites-only VCF, e.g. from vcfs_combine_batches
    training_intervals  // channel: [ path ] -- static resource, e.g. params.depth_training_bed (plain BED, not bgzipped -- GATK wants its own .idx, not a tabix .tbi; confirmed empirically, see README)
    median_coverage     // channel: [ meta, median_cov_bed ]
    rd_file             // channel: [ meta, rd_txt_gz, rd_txt_gz_tbi ] -- panel-wide, from vcfs_merge_read_counts
    pe_file             // channel: [ meta, pe_txt_gz, pe_txt_gz_tbi ] -- panel-wide, from vcfs_merge_evidence
    sr_file             // channel: [ meta, sr_txt_gz, sr_txt_gz_tbi ] -- panel-wide
    reference_dict      // channel: [ meta, dict ]
    ploidy_table        // channel: [ meta, ploidy_tsv ]
    depth_exclusion_intervals  // channel: [ bed, bed_tbi ]
    pesr_exclusion_intervals   // channel: [ bed, bed_tbi ]
    rf_cutoffs          // channel: [ meta, cutoffs_tsv ] -- from FilterBatchSites/AdjudicateSV
    contig_list         // channel: [ meta, contig_list ] -- plain-text, one contig per line
    cohort_name         // val: String, used as the output prefix
    extract_format_table_script  // val: absolute path to bin/extract_format_table.py

    main:
    GATK_TRAINSVGENOTYPING(
        vcf,
        training_intervals,
        median_coverage.map { meta, m -> m }.collect(),
        rd_file,
        pe_file,
        sr_file,
        reference_dict.map { meta, d -> d }.collect(),
        ploidy_table,
        depth_exclusion_intervals,
        pesr_exclusion_intervals,
        rf_cutoffs
    )

    contigs = contig_list
        .flatMap { meta, c -> file(c).readLines().findAll { it } }

    // Broadcast every "single panel-wide value" input across the
    // per-contig scatter below -- same wrap-collect-unwrap pattern as
    // vcfs_generate_batch_metrics (see that subworkflow's own note on why
    // .collect() alone flattens a channel of multi-element tuples).
    rd_value = rd_file.map { meta, f, tbi -> [ [ f, tbi ] ] }.collect().map { it[0] }
    pe_value = pe_file.map { meta, f, tbi -> [ [ f, tbi ] ] }.collect().map { it[0] }
    sr_value = sr_file.map { meta, f, tbi -> [ [ f, tbi ] ] }.collect().map { it[0] }
    dict_value = reference_dict.map { meta, d -> d }.collect()
    ploidy_value = ploidy_table.map { meta, t -> [ [ meta, t ] ] }.collect().map { it[0] }
    rd_table_value = GATK_TRAINSVGENOTYPING.out.rd_table.map { meta, t -> [ [ meta, t ] ] }.collect().map { it[0] }
    pe_table_value = GATK_TRAINSVGENOTYPING.out.pe_table.map { meta, t -> [ [ meta, t ] ] }.collect().map { it[0] }
    sr_table_value = GATK_TRAINSVGENOTYPING.out.sr_table.map { meta, t -> [ [ meta, t ] ] }.collect().map { it[0] }
    depth_excl_value = depth_exclusion_intervals.map { bed, tbi -> [ [ bed, tbi ] ] }.collect().map { it[0] }
    pesr_excl_value = pesr_exclusion_intervals.map { bed, tbi -> [ [ bed, tbi ] ] }.collect().map { it[0] }

    // meta.id = the contig name itself -- see top-of-file note on why
    // (joining the three independent PrintSVEvidence scatters back
    // together by a shared, meaningful key, not scatter/gather order).
    // .combine() flattens rd_value's wrapped [[f, tbi]] against the
    // scattered contig side too (same behavior as everywhere else in
    // this pipeline that combines a collected list -- confirmed here by
    // an actual arity mismatch at runtime, not assumed), so the
    // resulting closure takes 3 args (contig, f, tbi), not 2.
    per_contig_rd = contigs.combine(rd_value).map { contig, f, tbi -> [ [id: contig], f, tbi ] }
    per_contig_pe = contigs.combine(pe_value).map { contig, f, tbi -> [ [id: contig], f, tbi ] }
    per_contig_sr = contigs.combine(sr_value).map { contig, f, tbi -> [ [id: contig], f, tbi ] }

    GATK_PRINT_SV_EVIDENCE_CONTIG_RD(per_contig_rd, dict_value, contigs)
    GATK_PRINT_SV_EVIDENCE_CONTIG_PE(per_contig_pe, dict_value, contigs)
    GATK_PRINT_SV_EVIDENCE_CONTIG_SR(per_contig_sr, dict_value, contigs)

    TABIX_CONTIG_RD(GATK_PRINT_SV_EVIDENCE_CONTIG_RD.out.evidence)
    TABIX_CONTIG_PE(GATK_PRINT_SV_EVIDENCE_CONTIG_PE.out.evidence)
    TABIX_CONTIG_SR(GATK_PRINT_SV_EVIDENCE_CONTIG_SR.out.evidence)

    // Join the three per-contig evidence channels back together by their
    // shared contig-name key (meta.id), then rebuild the
    // "${cohort_name}.genotype_batch.${contig}" naming GenotypeSVs' own
    // output filename should carry (matching the WDL's own
    // output_prefix) and the -L <contig> value GenotypeSVs itself needs.
    // The sites VCF is the same single value for every contig -- broadcast
    // it in via .combine(), same pattern as every other "single value
    // across a scatter" case in this subworkflow.
    vcf_value = vcf.map { meta, v, tbi -> [ [ v, tbi ] ] }.collect().map { it[0] }

    rd_pe_sr_by_contig = TABIX_CONTIG_RD.out.evidence
        .map { meta, f, tbi -> [ meta.id, f, tbi ] }
        .join(TABIX_CONTIG_PE.out.evidence.map { meta, f, tbi -> [ meta.id, f, tbi ] })
        .join(TABIX_CONTIG_SR.out.evidence.map { meta, f, tbi -> [ meta.id, f, tbi ] })

    // Same .combine() flattening as per_contig_rd/pe/sr above -- vcf_value
    // ([v, tbi]) flattens against rd_pe_sr_by_contig's own 7 positional
    // elements, giving 9 total (not 8, "..., vcf_pair" as one element).
    genotype_svs_input = rd_pe_sr_by_contig
        .combine(vcf_value)
        .map { contig, rd, rd_tbi, pe, pe_tbi, sr, sr_tbi, v, vcf_tbi ->
            [
                [ [id: "${cohort_name}.genotype_batch.${contig}"], v, vcf_tbi ],
                contig,
                [ [id: contig], rd, rd_tbi ],
                [ [id: contig], pe, pe_tbi ],
                [ [id: contig], sr, sr_tbi ]
            ]
        }

    GATK_GENOTYPESVS(
        genotype_svs_input.map { v, contig, rd, pe, sr -> v },
        genotype_svs_input.map { v, contig, rd, pe, sr -> contig },
        median_coverage.map { meta, m -> m }.collect(),
        genotype_svs_input.map { v, contig, rd, pe, sr -> rd },
        genotype_svs_input.map { v, contig, rd, pe, sr -> pe },
        genotype_svs_input.map { v, contig, rd, pe, sr -> sr },
        dict_value,
        ploidy_value,
        depth_excl_value,
        pesr_excl_value,
        rd_table_value,
        pe_table_value,
        sr_table_value
    )

    concat_input = GATK_GENOTYPESVS.out.vcf
        .map { meta, v, tbi -> [ v, tbi ] }
        .collect(flat: false)
        .map { pairs -> pairs.sort { it[0].name } }
        .map { pairs -> [ [id: "${cohort_name}.genotype_batch"], pairs.collect { it[0] }, pairs.collect { it[1] } ] }

    CONCAT_GENOTYPE_BATCH_SHARDS(concat_input)

    SEPARATE_DEPTH_PESR(CONCAT_GENOTYPE_BATCH_SHARDS.out.vcf)

    GENERATE_REGENO_COVERAGE_MEDIANS(
        SEPARATE_DEPTH_PESR.out.depth_vcf.map { meta, v, tbi -> [ [id: "${cohort_name}.regeno_coverage_medians"], v, tbi ] },
        file(extract_format_table_script)
    )

    emit:
    genotyped_depth_vcf = SEPARATE_DEPTH_PESR.out.depth_vcf
    genotyped_pesr_vcf = SEPARATE_DEPTH_PESR.out.pesr_vcf
    genotyping_rd_table = GATK_TRAINSVGENOTYPING.out.rd_table
    genotyping_pe_table = GATK_TRAINSVGENOTYPING.out.pe_table
    genotyping_sr_table = GATK_TRAINSVGENOTYPING.out.sr_table
    regeno_coverage_medians = GENERATE_REGENO_COVERAGE_MEDIANS.out.regeno_coverage_medians
}
