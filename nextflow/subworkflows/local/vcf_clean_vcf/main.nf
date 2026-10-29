//
// Final cleanup of the cohort VCF: sex-chromosome genotypes, copy-number
// consistency between overlapping CNVs, multiallelic CNVs, redundant and
// fragmented CNVs, mobile-element deletions, low-precision filters, final
// IDs. Mirrors GATK-SV's CleanVcf (wdl/CleanVcf.wdl, with
// wdl/CleanVcfChromosome.wdl), the last step of MakeCohortVcf, on
// GenotypeComplexVariants' output.
//
// Per contig: CLEANVCF_PREPROCESS -> CLEANVCF_REVISE -> CLEANVCF_POSTPROCESS
// (upstream's ~15 tasks in 3, grouped by container; contigs unsharded,
// upstream: 5,000 records per shard). Then a genome-wide concat. Output is
// GATK format, as upstream.
//
// Not implemented: the outlier-sample list (no FilterBatchSamples yet) and
// MakeCohortVcfMetrics (off by default upstream).
//

include { PED_SEXES; pedSexes  } from '../../../modules/local/ped_sexes/main'
include { CLEANVCF_PREPROCESS  } from '../../../modules/local/clean_vcf/preprocess/main'
include { CLEANVCF_REVISE      } from '../../../modules/local/clean_vcf/revise/main'
include { CLEANVCF_POSTPROCESS } from '../../../modules/local/clean_vcf/postprocess/main'
include { CONCAT_VCFS as CONCAT_CLEANED_VCFS } from '../../../modules/local/concat_vcfs/main'

workflow VCF_CLEAN_VCF {
    take:
    vcf                 // channel: [ meta, vcf, vcf_tbi ] -- GenotypeComplexVariants' complex_genotype VCF
    bothside_pass       // channel: [ meta, txt ] -- ResolveComplexVariants' updated BOTHSIDES_SUPPORT list
    background_fail     // channel: [ meta, txt ] -- ResolveComplexVariants' updated HIGH_SR_BACKGROUND list
    ped                 // channel: [ meta, ped ]
    ploidy_table        // channel: [ meta, tsv ] -- retain_female_chr_y=false, as upstream's FormatVcfForGatk
    contig_list         // channel: [ meta, contig_list ]
    mei_references      // val: [ line1_bed_gz, hervk_bed_gz, intron_tsv_gz ]
    chr_x               // val: String
    chr_y               // val: String
    cohort_name         // val: String
    preprocess_script   // val: absolute path to bin/cleanvcf_preprocess.py
    retro_del_script    // val: absolute path to bin/add_retro_del_filters.py

    main:
    contigs = contig_list.flatMap { meta, c -> file(c).readLines().findAll { it } }

    PED_SEXES(ped.map { meta, f -> [ meta, pedSexes(f) ] })
    sexes_value      = PED_SEXES.out.ped.map { meta, p -> p }.collect()
    ploidy_value     = ploidy_table.map { meta, t -> t }.collect()
    bothside_value   = bothside_pass.map { meta, f -> f }.collect()
    background_value = background_fail.map { meta, f -> f }.collect()

    CLEANVCF_PREPROCESS(
        vcf.combine(contigs)
            .map { meta, v, tbi, contig -> [ [id: "${cohort_name}.${contig}.preprocess", contig: contig, cohort: cohort_name], v, tbi ] },
        ploidy_value, bothside_value, background_value, sexes_value, chr_x, chr_y,
        file(preprocess_script)
    )

    CLEANVCF_REVISE(
        CLEANVCF_PREPROCESS.out.vcf
            .map { meta, v, tbi -> [ meta + [id: "${cohort_name}.${meta.contig}.revise"], v, tbi ] }
    )

    CLEANVCF_POSTPROCESS(
        CLEANVCF_REVISE.out.vcf
            .map { meta, v, tbi -> [ meta + [id: "${cohort_name}.${meta.contig}.cleaned"], v, tbi ] },
        ploidy_value, bothside_value, background_value, sexes_value, mei_references,
        file(retro_del_script)
    )

    CONCAT_CLEANED_VCFS(
        CLEANVCF_POSTPROCESS.out.vcf
            .map { meta, v, tbi -> [ v, tbi ] }
            .collect(flat: false)
            .map { pairs -> pairs.sort { it[0].name } }
            .map { pairs -> [ [id: "${cohort_name}.cleaned"], pairs.collect { it[0] }, pairs.collect { it[1] } ] }
    )

    emit:
    vcf = CONCAT_CLEANED_VCFS.out.vcf   // channel: [ meta, vcf, vcf_tbi ]
}
