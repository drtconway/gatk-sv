#!/usr/bin/env nextflow

//
// Standalone entry point: sample sheet -> genotyped depth/PESR VCFs --
// the actual goal of the genotyping prerequisite chain (see
// ../README.md#path-to-genotyping-the-real-prerequisite-chain).
//
// Run directly, e.g.:
//   nextflow run workflows/genotype_batch.nf -profile apptainer \
//     --input samplesheet.tsv --fasta ref.fa --fasta_fai ref.fa.fai \
//     --primary_contigs_list primary_contigs.list \
//     --primary_contigs_fai contigs.fai --ped cohort.ped \
//     --pesr_exclude_intervals exclude.bed.gz \
//     --pesr_exclude_intervals_tbi exclude.bed.gz.tbi \
//     --manta_region_bed manta_region.bed.gz \
//     --manta_region_bed_tbi manta_region.bed.gz.tbi \
//     --reference_dict ref.dict \
//     --sd_locs_vcf dbsnp.vcf --sd_locs_vcf_idx dbsnp.vcf.idx \
//     --preprocessed_intervals preprocessed_intervals.interval_list \
//     --segdups segdups.bed.gz --segdups_tbi segdups.bed.gz.tbi \
//     --rmsk rmsk.bed.gz --rmsk_tbi rmsk.bed.gz.tbi \
//     --depth_training_bed train.bed \
//     --bin_exclude bin_exclude.bed.gz --bin_exclude_tbi bin_exclude.bed.gz.tbi \
//     --pesr_exclude_list pesr_exclude.bed.gz --pesr_exclude_list_tbi pesr_exclude.bed.gz.tbi \
//     --outdir results
//
// Composes generate_batch_metrics.nf's own full chain (stage-1
// clustering, panel-wide evidence, GenerateBatchMetrics) with
// AdjudicateSV (rf_cutoffs), the new filter_batch_sites (RF-score sites
// filtering + cross-caller merge), and genotype_batch itself -- same
// "one entry point runs the full chain" pattern as every other
// workflows/*.nf entry point. See
// README.md#path-to-genotyping-the-real-prerequisite-chain for why this
// is the last of three sequential pieces (GenerateBatchMetrics ->
// FilterBatchSites/AdjudicateSV -> GenotypeBatch).
//
// depth_training_bed is deliberately a plain (uncompressed) BED, not
// bgzipped/tabix-indexed like most other static resources here --
// confirmed empirically (GATK's TrainSVGenotyping wants its own .idx for
// a bgzipped input, not a tabix .tbi; see subworkflows/local/
// genotype_batch's own note), matching the real upstream resource shape
// (inputs/values/resources_hg38.json's depth_training_bed is also plain
// .bed, not .bed.gz).
//

include { UTILS_INPUT_CHANNELS        } from '../subworkflows/local/utils_input_channels/main'
include { BAM_CALL_MANTA              } from '../subworkflows/local/bam_call_manta/main'
include { BAM_CALL_WHAM               } from '../subworkflows/local/bam_call_wham/main'
include { VCFS_CLUSTER_SVCLUSTER as VCFS_CLUSTER_SVCLUSTER_MANTA } from '../subworkflows/local/vcfs_cluster_svcluster/main'
include { VCFS_CLUSTER_SVCLUSTER as VCFS_CLUSTER_SVCLUSTER_WHAM  } from '../subworkflows/local/vcfs_cluster_svcluster/main'
include { BAM_COLLECT_EVIDENCE        } from '../subworkflows/local/bam_collect_evidence/main'
include { VCFS_MERGE_EVIDENCE         } from '../subworkflows/local/vcfs_merge_evidence/main'
include { VCFS_MERGE_READ_COUNTS      } from '../subworkflows/local/vcfs_merge_read_counts/main'
include { MEDIAN_COVERAGE             } from '../modules/local/median_coverage/main'
include { VCFS_GENERATE_BATCH_METRICS } from '../subworkflows/local/vcfs_generate_batch_metrics/main'
include { ADJUDICATE_SV               } from '../modules/local/adjudicate_sv/main'
include { FILTER_BATCH_SITES          } from '../subworkflows/local/filter_batch_sites/main'
include { GENOTYPE_BATCH              } from '../subworkflows/local/genotype_batch/main'

workflow {
    main:
    def pipelineRoot = "${projectDir}/.."

    [
        'primary_contigs_list', 'primary_contigs_fai', 'ped',
        'pesr_exclude_intervals', 'pesr_exclude_intervals_tbi',
        'reference_dict', 'manta_region_bed', 'manta_region_bed_tbi',
        'sd_locs_vcf', 'sd_locs_vcf_idx', 'preprocessed_intervals',
        'segdups', 'segdups_tbi', 'rmsk', 'rmsk_tbi',
        'depth_training_bed', 'bin_exclude', 'bin_exclude_tbi',
        'pesr_exclude_list', 'pesr_exclude_list_tbi'
    ].each { p ->
        if (!params[p]) {
            error "genotype_batch.nf requires --${p}"
        }
    }

    UTILS_INPUT_CHANNELS(params.input, pipelineRoot, params.fasta, params.fasta_fai)
    samples   = UTILS_INPUT_CHANNELS.out.samples
    fasta     = UTILS_INPUT_CHANNELS.out.fasta
    fasta_fai = UTILS_INPUT_CHANNELS.out.fasta_fai

    contig_list = channel.of([ [id: 'contigs'], file(params.primary_contigs_list) ])
    contigs_fai = channel.of([ [id: 'contigs'], file(params.primary_contigs_fai) ])
    manta_region_bed = channel.of([
        [id: 'manta_region'],
        file(params.manta_region_bed),
        file(params.manta_region_bed_tbi)
    ])
    ped = channel.of([ [id: 'ped'], file(params.ped) ])
    exclude_intervals = channel.of([
        [id: 'exclude'],
        file(params.pesr_exclude_intervals),
        file(params.pesr_exclude_intervals_tbi)
    ])
    dict = channel.of([ [id: 'reference'], file(params.reference_dict) ])
    sd_locs_vcf = channel.of([
        [id: 'sd_locs'],
        file(params.sd_locs_vcf),
        file(params.sd_locs_vcf_idx)
    ])
    preprocessed_intervals = channel.of([ [id: 'intervals'], file(params.preprocessed_intervals) ])
    segdups = channel.of([ [id: 'segdups'], file(params.segdups), file(params.segdups_tbi) ])
    rmsk = channel.of([ [id: 'rmsk'], file(params.rmsk), file(params.rmsk_tbi) ])

    // --- Per-caller call + cluster (stage 1), same as combine_batches.nf ---

    BAM_CALL_MANTA(samples, fasta, fasta_fai, manta_region_bed, contigs_fai, params.min_svsize)
    VCFS_CLUSTER_SVCLUSTER_MANTA(
        BAM_CALL_MANTA.out.vcf,
        ped, contig_list, exclude_intervals, fasta, fasta_fai, dict,
        'manta', params.min_svsize,
        "${pipelineRoot}/bin/ploidy_table_from_ped.py",
        "${pipelineRoot}/bin/format_svtk_vcf_for_gatk.py",
        "${pipelineRoot}/bin/format_gatk_vcf_for_svtk.py"
    )

    BAM_CALL_WHAM(samples, fasta, fasta_fai, contigs_fai, params.min_svsize)
    VCFS_CLUSTER_SVCLUSTER_WHAM(
        BAM_CALL_WHAM.out.vcf,
        ped, contig_list, exclude_intervals, fasta, fasta_fai, dict,
        'wham', params.min_svsize,
        "${pipelineRoot}/bin/ploidy_table_from_ped.py",
        "${pipelineRoot}/bin/format_svtk_vcf_for_gatk.py",
        "${pipelineRoot}/bin/format_gatk_vcf_for_svtk.py"
    )

    caller_vcfs = VCFS_CLUSTER_SVCLUSTER_MANTA.out.clustered_vcf
        .mix(VCFS_CLUSTER_SVCLUSTER_WHAM.out.clustered_vcf)

    // --- Panel-wide evidence, same as merge_evidence.nf ---

    BAM_COLLECT_EVIDENCE(
        samples,
        fasta,
        fasta_fai,
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

    MEDIAN_COVERAGE(
        VCFS_MERGE_READ_COUNTS.out.merged_counts.map { meta, rd, tbi -> [ meta, rd ] },
        file("${pipelineRoot}/bin/medianCoverage.R")
    )

    // --- Per-variant evidence metrics ---

    VCFS_GENERATE_BATCH_METRICS(
        caller_vcfs,
        ped,
        contig_list,
        dict,
        VCFS_MERGE_EVIDENCE.out.merged_pe,
        VCFS_MERGE_EVIDENCE.out.merged_sr,
        VCFS_MERGE_EVIDENCE.out.merged_baf,
        VCFS_MERGE_READ_COUNTS.out.merged_counts,
        MEDIAN_COVERAGE.out.median_coverage,
        segdups,
        rmsk,
        params.chr_x ?: 'chrX',
        params.chr_y ?: 'chrY',
        'cohort',
        (params.records_per_shard ?: 10000) as Integer,
        "${pipelineRoot}/bin/ploidy_table_from_ped.py",
        "${pipelineRoot}/bin/format_svtk_vcf_for_gatk.py",
        "${pipelineRoot}/bin/aggregate.py"
    )

    // --- Random-forest cutoff derivation ---

    ADJUDICATE_SV(
        VCFS_GENERATE_BATCH_METRICS.out.metrics,
        file("${pipelineRoot}/bin/adjudicate.py"),
        file("${pipelineRoot}/bin/adjudicate_sv.py"),
        file("${pipelineRoot}/bin/random_forest.py"),
        file("${pipelineRoot}/bin/labelers.py")
    )

    // --- RF-score sites filtering + cross-caller merge ---

    FILTER_BATCH_SITES(
        caller_vcfs,
        VCFS_GENERATE_BATCH_METRICS.out.metrics,
        ADJUDICATE_SV.out.scores,
        ADJUDICATE_SV.out.cutoffs,
        'cohort',
        "${pipelineRoot}/bin/rewrite_SR_coords.py",
        "${pipelineRoot}/bin/annotate_RF_evidence.py"
    )

    // --- Genotyping ---

    depth_exclusion_intervals = channel.value([
        file(params.bin_exclude), file(params.bin_exclude_tbi)
    ])
    pesr_exclusion_intervals = channel.value([
        file(params.pesr_exclude_list), file(params.pesr_exclude_list_tbi)
    ])
    training_intervals = channel.value(file(params.depth_training_bed))

    GENOTYPE_BATCH(
        FILTER_BATCH_SITES.out.merged_vcf,
        training_intervals,
        MEDIAN_COVERAGE.out.median_coverage,
        VCFS_MERGE_READ_COUNTS.out.merged_counts,
        VCFS_MERGE_EVIDENCE.out.merged_pe,
        VCFS_MERGE_EVIDENCE.out.merged_sr,
        dict,
        VCFS_GENERATE_BATCH_METRICS.out.ploidy_table.map { list -> [ [id: 'ploidy'], list[0] ] },
        depth_exclusion_intervals,
        pesr_exclusion_intervals,
        ADJUDICATE_SV.out.cutoffs,
        contig_list,
        'cohort',
        "${pipelineRoot}/bin/extract_format_table.py"
    )

    GENOTYPE_BATCH.out.genotyped_depth_vcf.view { meta, v, i -> "Genotyped depth VCF: ${v}" }
    GENOTYPE_BATCH.out.genotyped_pesr_vcf.view { meta, v, i -> "Genotyped PESR VCF: ${v}" }
}
