//
// Depth-genotype complex SVs: test each CPX record's (and each unresolved
// inversion single-ender's) CNV intervals with RdTest across all samples,
// then confirm, reclassify or drop the record. Mirrors GATK-SV's
// GenotypeComplexVariants (wdl/GenotypeComplexVariants.wdl, with
// ScatterCpxGenotyping/GenotypeCpxCnvs/GenotypeCpxCnvsPerBatch), the third
// step of MakeCohortVcf, on ResolveComplexVariants' output.
//
// Per contig, as upstream:
//   1. CPX_CNV_INTERVALS: the contig's records and their CNV intervals
//   2. SPLIT_BED_BY_SIZE: chunk the intervals
//   3. RDTEST_GENOTYPE per chunk
//   4. PARSE_CPX_GENOTYPES: merge the chunks' genotypes and apply them
// Then a genome-wide concat.
//
// Single batch, so upstream's per-batch scatter and PED subsetting go
// away, and the contig isn't sharded further (upstream: 50,000 records per
// shard; a contig here has a few thousand).
//

include { PED_SEXES; pedSexes  } from '../../../modules/local/ped_sexes/main'
include { CPX_CNV_INTERVALS    } from '../../../modules/local/genotype_complex/cpx_cnv_intervals/main'
include { SPLIT_BED_BY_SIZE    } from '../../../modules/local/genotype_complex/split_bed_by_size/main'
include { RDTEST_GENOTYPE      } from '../../../modules/local/genotype_complex/rdtest_genotype/main'
include { PARSE_CPX_GENOTYPES  } from '../../../modules/local/genotype_complex/parse_cpx_genotypes/main'
include { CONCAT_VCFS as CONCAT_COMPLEX_GENOTYPE } from '../../../modules/local/concat_vcfs/main'

workflow VCF_GENOTYPE_COMPLEX_VARIANTS {
    take:
    vcf                 // channel: [ meta, vcf, vcf_tbi ] -- ResolveComplexVariants' complex_resolve VCF
    ped                 // channel: [ meta, ped ]
    coverage_file       // channel: [ meta, rd_txt_gz, tbi ] -- panel-wide read-count matrix (vcfs_merge_read_counts)
    median_coverage     // channel: [ meta, median_cov_bed ]
    rd_table            // channel: [ meta, tsv ] -- TrainSVGenotyping's rd_table
    bin_exclude         // val: [ bed_gz, tbi ]
    contig_list         // channel: [ meta, contig_list ]
    reference_dict      // channel: [ meta, dict ]
    cohort_name         // val: String

    main:
    contigs = contig_list.flatMap { meta, c -> file(c).readLines().findAll { it } }

    PED_SEXES(ped.map { meta, f -> [ meta, pedSexes(f) ] })
    sexes_value    = PED_SEXES.out.ped.map { meta, p -> p }.collect()
    coverage_value = coverage_file.map { meta, f, tbi -> [ [ f, tbi ] ] }.collect().map { it[0] }
    median_value   = median_coverage.map { meta, m -> m }.collect()
    rd_table_value = rd_table.map { meta, t -> t }.collect()
    dict_value     = reference_dict.map { meta, d -> d }.collect()

    // --- 1. per contig ---
    CPX_CNV_INTERVALS(
        vcf.combine(contigs)
            .map { meta, v, tbi, contig -> [ [id: "${cohort_name}.${contig}", contig: contig], v, tbi ] }
    )

    // --- 2-3. chunk and genotype; upstream's split sizes and RdTest bins ---
    SPLIT_BED_BY_SIZE(CPX_CNV_INTERVALS.out.intervals, median_value, 2500, 250)

    chunks = SPLIT_BED_BY_SIZE.out.beds
        .flatMap { meta, beds ->
            (beds instanceof List ? beds : [ beds ]).collect { b -> [ meta + [id: "${meta.id}.${b.name}", key: meta.id], b ] }
        }

    RDTEST_GENOTYPE(chunks, coverage_value, median_value, sexes_value, bin_exclude, rd_table_value, dict_value, 100000)

    // --- 4. per contig, chunks sorted for a stable hash ---
    genotypes = RDTEST_GENOTYPE.out.genotypes
        .map { meta, g -> [ meta.key, g ] }
        .groupTuple()
        .map { key, gs -> [ key, gs.sort { it.name } ] }

    parse_input = CPX_CNV_INTERVALS.out.vcf.map { meta, v, tbi -> [ meta.id, meta, v, tbi ] }
        .join(CPX_CNV_INTERVALS.out.intervals.map { meta, b -> [ meta.id, b ] })
        .join(genotypes)
        .map { key, meta, v, tbi, intervals, gs -> [ meta + [id: "${meta.id}.complex_genotype"], v, tbi, intervals, gs ] }

    PARSE_CPX_GENOTYPES(parse_input, sexes_value)

    // --- genome-wide ---
    CONCAT_COMPLEX_GENOTYPE(
        PARSE_CPX_GENOTYPES.out.vcf
            .map { meta, v, tbi -> [ v, tbi ] }
            .collect(flat: false)
            .map { pairs -> pairs.sort { it[0].name } }
            .map { pairs -> [ [id: "${cohort_name}.complex_genotype"], pairs.collect { it[0] }, pairs.collect { it[1] } ] }
    )

    emit:
    vcf                    = CONCAT_COMPLEX_GENOTYPE.out.vcf                     // channel: [ meta, vcf, vcf_tbi ]
    reclassification_table = PARSE_CPX_GENOTYPES.out.reclassification_table      // channel: [ meta, txt_gz ] -- per contig
}
