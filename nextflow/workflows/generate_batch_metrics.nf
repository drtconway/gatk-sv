#!/usr/bin/env nextflow

//
// Standalone entry point: sample sheet -> per-variant evidence metrics
// (evidence_metrics.tsv), the input FilterBatchSites/AdjudicateSV needs
// to derive rf_cutoffs for genotyping -- see
// ../README.md#path-to-genotyping-the-real-prerequisite-chain.
//
// Run directly, e.g.:
//   nextflow run workflows/generate_batch_metrics.nf -profile apptainer \
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
//     --outdir results
//
// Composes bam_call_manta/bam_call_wham + vcfs_cluster_svcluster (stage 1
// per-caller clustering, same as combine_batches.nf's own first half) with
// bam_collect_evidence + vcfs_merge_evidence + vcfs_merge_read_counts
// (panel-wide evidence, same as merge_evidence.nf) and the new
// vcfs_generate_batch_metrics -- same "one entry point runs the full
// chain" pattern as combine_batches.nf/merge_evidence.nf. Deliberately
// does NOT run vcfs_combine_batches: GenerateBatchMetrics takes each
// caller's stage-1 clustered VCF directly, not CombineBatches' later
// cross-caller merged output -- confirmed against
// wdl/GATKSVPipelinePhase1.wdl's own wiring (GenerateBatchMetrics.
// {depth,manta,wham,...}_vcf = ClusterBatch.clustered_*_vcf), not assumed
// from the workflow name alone. See
// subworkflows/local/vcfs_generate_batch_metrics's own top-of-file note.
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

workflow {
    main:
    def pipelineRoot = "${projectDir}/.."

    [
        'primary_contigs_list', 'primary_contigs_fai', 'ped',
        'pesr_exclude_intervals', 'pesr_exclude_intervals_tbi',
        'reference_dict', 'manta_region_bed', 'manta_region_bed_tbi',
        'sd_locs_vcf', 'sd_locs_vcf_idx', 'preprocessed_intervals',
        'segdups', 'segdups_tbi', 'rmsk', 'rmsk_tbi'
    ].each { p ->
        if (!params[p]) {
            error "generate_batch_metrics.nf requires --${p}"
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

    VCFS_GENERATE_BATCH_METRICS.out.metrics.view { meta, tsv -> "Evidence metrics: ${tsv}" }
}
