//
// Cluster one caller's SV calls across all samples in a run into
// representative sites, mirroring GATK-SV's ClusterPESR workflow
// (wdl/PESRClustering.wdl) for PE/SR-based callers (Manta, Wham,
// Scramble). This is stage 1 of harmonisation -- per-caller, cross-sample
// clustering -- as distinct from stage 2, cross-caller merging
// (CombineBatches/MergeBatchSites), which is not implemented yet. See
// README.md#how-gatk-svs-harmonisation-actually-works.
//
// Steps, matching ClusterPESR:
//   1. Ploidy table from the pedigree file (CreatePloidyTableFromPed)
//   2. Per-sample: format conversion + interval/size filtering
//      (PreparePESRVcfs)
//   3. SVCluster across all samples' prepared VCFs in one call
//      (not scattered by contig -- see note below)
//   4. Post-clustering interval exclusion (ExcludeIntervalsByEndpoints)
//   5. Format back to svtk (GatkToSvtkVcf)
//
// Deliberately not implemented: GATK-SV scatters step 3 by contig then
// concatenates, for parallelism at cohort scale. We run one SVCluster
// call across the whole genome instead -- revisit if runtime on real
// data warrants it (same simplification already made for Wham's region
// scatter, see subworkflows/local/bam_call_wham).
//

include { GENOME_FILE              } from '../../../modules/local/genome_file/main'
include { PLOIDY_TABLE_FROM_PED    } from '../../../modules/local/ploidy_table_from_ped/main'
include { FORMAT_SVTK_VCF_FOR_GATK } from '../../../modules/local/format_svtk_vcf_for_gatk/main'
include { VCF_ENDS_BED as VCF_ENDS_BED_PRE  } from '../../../modules/local/vcf_ends_bed/main'
include { VCF_ENDS_BED as VCF_ENDS_BED_POST } from '../../../modules/local/vcf_ends_bed/main'
include { EXCLUDED_VARIANT_IDS as EXCLUDED_VARIANT_IDS_PRE  } from '../../../modules/local/excluded_variant_ids/main'
include { EXCLUDED_VARIANT_IDS as EXCLUDED_VARIANT_IDS_POST } from '../../../modules/local/excluded_variant_ids/main'
include { BEDTOOLS_INTERSECT as BEDTOOLS_INTERSECT_PRE  } from '../../../modules/nf-core/bedtools/intersect/main'
include { BEDTOOLS_INTERSECT as BEDTOOLS_INTERSECT_POST } from '../../../modules/nf-core/bedtools/intersect/main'
include { EXCLUDE_VARIANTS_BY_ID as EXCLUDE_VARIANTS_BY_ID_PRE  } from '../../../modules/local/exclude_variants_by_id/main'
include { EXCLUDE_VARIANTS_BY_ID as EXCLUDE_VARIANTS_BY_ID_POST } from '../../../modules/local/exclude_variants_by_id/main'
include { GATK4_SVCLUSTER          } from '../../../modules/nf-core/gatk4/svcluster/main'
include { FORMAT_GATK_VCF_FOR_SVTK } from '../../../modules/local/format_gatk_vcf_for_svtk/main'

workflow VCFS_CLUSTER_SVCLUSTER {
    take:
    vcfs               // channel: [ meta, vcf ] -- one per sample, same caller (meta.id = sample_id)
    ped                // channel: [ meta, ped ] -- one pedigree file for the whole run
    contig_list        // channel: [ meta, contig_list ]
    exclude_intervals  // channel: [ meta, bed, bed_tbi ] -- e.g. params.pesr_exclude_intervals
    fasta              // channel: [ meta, fasta ]
    fasta_fai          // channel: [ meta, fasta_fai ]
    dict               // channel: [ meta, dict ] -- sequence dictionary, required by GATK4_SVCLUSTER
    caller              // val: caller name, e.g. 'manta' -- used as --source for the final format conversion
    min_size            // val: Int, minimum SV size to retain (GATK-SV default: 50) -- not yet applied, see note below

    main:
    GENOME_FILE(fasta_fai)
    genome_file = GENOME_FILE.out.genome_file.map { meta, g -> g }.collect()

    PLOIDY_TABLE_FROM_PED(ped, contig_list.map { meta, c -> c }.collect())
    ploidy_table = PLOIDY_TABLE_FROM_PED.out.ploidy_table.map { meta, t -> t }.collect()

    // --- Step 2: per-sample prepare (format conversion + interval filter) ---
    // NB: GATK-SV's min_size filtering (SVLEN>=min_size, via bcftools view
    // -i alongside the ID exclusion) isn't applied yet -- only interval
    // exclusion is. Revisit once this subworkflow's overall shape is
    // validated; min_size is threaded through as a param already.

    FORMAT_SVTK_VCF_FOR_GATK(vcfs, ploidy_table.map { t -> [ [id: 'ploidy'], t ] })

    VCF_ENDS_BED_PRE(FORMAT_SVTK_VCF_FOR_GATK.out.vcf)

    intersect_input_pre = VCF_ENDS_BED_PRE.out.bed
        .combine(exclude_intervals.map { meta, bed, tbi -> bed })
        .map { meta, ends_bed, exclude_bed -> [ meta, ends_bed, exclude_bed ] }
    BEDTOOLS_INTERSECT_PRE(intersect_input_pre, genome_file.map { g -> [ [id: 'genome'], g ] })

    EXCLUDED_VARIANT_IDS_PRE(BEDTOOLS_INTERSECT_PRE.out.intersect)

    EXCLUDE_VARIANTS_BY_ID_PRE(
        FORMAT_SVTK_VCF_FOR_GATK.out.vcf.join(EXCLUDED_VARIANT_IDS_PRE.out.excluded_ids)
    )
    prepared_vcfs = EXCLUDE_VARIANTS_BY_ID_PRE.out.vcf

    // --- Step 3: SVCluster across all samples in one call ---

    cluster_input = prepared_vcfs
        .map { meta, vcf, tbi -> [ vcf, tbi ] }
        .collect(flat: false)
        .map { pairs -> [ [id: caller], pairs.collect { it[0] }, pairs.collect { it[1] } ] }

    GATK4_SVCLUSTER(
        cluster_input,
        ploidy_table,
        fasta.map { meta, f -> f }.collect(),
        fasta_fai.map { meta, f -> f }.collect(),
        dict.map { meta, d -> d }.collect()
    )

    // --- Step 4: post-clustering interval exclusion ---

    VCF_ENDS_BED_POST(GATK4_SVCLUSTER.out.clustered_vcf)

    intersect_input_post = VCF_ENDS_BED_POST.out.bed
        .combine(exclude_intervals.map { meta, bed, tbi -> bed })
        .map { meta, ends_bed, exclude_bed -> [ meta, ends_bed, exclude_bed ] }
    BEDTOOLS_INTERSECT_POST(intersect_input_post, genome_file.map { g -> [ [id: 'genome'], g ] })

    EXCLUDED_VARIANT_IDS_POST(BEDTOOLS_INTERSECT_POST.out.intersect)

    // GATK4_SVCLUSTER's clustered_vcf and clustered_vcf_index are
    // separate emit channels, not a single [meta,vcf,tbi] tuple, so join
    // them before feeding the [meta,vcf] shape EXCLUDE_VARIANTS_BY_ID
    // expects (its own vcf/tbi outputs aren't needed here downstream).
    clustered_vcf_only = GATK4_SVCLUSTER.out.clustered_vcf
    EXCLUDE_VARIANTS_BY_ID_POST(
        clustered_vcf_only.join(EXCLUDED_VARIANT_IDS_POST.out.excluded_ids)
    )

    // --- Step 5: format back to svtk ---

    FORMAT_GATK_VCF_FOR_SVTK(
        EXCLUDE_VARIANTS_BY_ID_POST.out.vcf.map { meta, vcf, tbi -> [ meta, vcf ] },
        contig_list.map { meta, c -> c }.collect(),
        caller
    )

    clustered_vcf = FORMAT_GATK_VCF_FOR_SVTK.out.vcf

    emit:
    clustered_vcf  // channel: [ meta, vcf, vcf_tbi ] -- meta.id == caller
}
