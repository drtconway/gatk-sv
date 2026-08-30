#!/usr/bin/env nextflow

//
// Standalone entry point: sample sheet -> per-caller clustering (Manta,
// Wham) -> one cross-caller merged site VCF (harmonisation stage 2 -- see
// subworkflows/local/vcfs_combine_batches).
//
// Run directly, e.g.:
//   nextflow run workflows/combine_batches.nf -profile apptainer \
//     --input samplesheet.tsv --fasta ref.fa --fasta_fai ref.fa.fai \
//     --primary_contigs_list primary_contigs.list \
//     --primary_contigs_fai contigs.fai --ped cohort.ped \
//     --pesr_exclude_intervals exclude.bed.gz \
//     --pesr_exclude_intervals_tbi exclude.bed.gz.tbi \
//     --manta_region_bed manta_region.bed.gz \
//     --manta_region_bed_tbi manta_region.bed.gz.tbi \
//     --reference_dict ref.dict \
//     --clustering_config_part1 clustering_config.part_one.tsv \
//     --stratification_config_part1 stratify_config.part_one.tsv \
//     --clustering_config_part2 clustering_config.part_two.tsv \
//     --stratification_config_part2 stratify_config.part_two.tsv \
//     --clustering_track_sr hg38.SimpRep.bed \
//     --clustering_track_sd hg38.SegDup.bed \
//     --clustering_track_rm hg38.RM.bed \
//     --outdir results
//
// One sample sheet in, one cross-caller merged site VCF out: composes the
// same call+cluster subworkflows cluster_manta.nf/cluster_wham.nf use
// standalone, then feeds both callers' clustered output into
// vcfs_combine_batches. See ../README.md#composability.
//

include { UTILS_INPUT_CHANNELS   } from '../subworkflows/local/utils_input_channels/main'
include { BAM_CALL_MANTA         } from '../subworkflows/local/bam_call_manta/main'
include { BAM_CALL_WHAM          } from '../subworkflows/local/bam_call_wham/main'
include { VCFS_CLUSTER_SVCLUSTER as VCFS_CLUSTER_SVCLUSTER_MANTA } from '../subworkflows/local/vcfs_cluster_svcluster/main'
include { VCFS_CLUSTER_SVCLUSTER as VCFS_CLUSTER_SVCLUSTER_WHAM  } from '../subworkflows/local/vcfs_cluster_svcluster/main'
include { VCFS_COMBINE_BATCHES   } from '../subworkflows/local/vcfs_combine_batches/main'

workflow {
    main:
    def pipelineRoot = "${projectDir}/.."

    [
        'primary_contigs_list', 'primary_contigs_fai', 'ped',
        'pesr_exclude_intervals', 'pesr_exclude_intervals_tbi',
        'reference_dict', 'manta_region_bed', 'manta_region_bed_tbi',
        'clustering_config_part1', 'stratification_config_part1',
        'clustering_config_part2', 'stratification_config_part2',
        'clustering_track_sr', 'clustering_track_sd', 'clustering_track_rm'
    ].each { p ->
        if (!params[p]) {
            error "combine_batches.nf requires --${p}"
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

    // --- Per-caller call + cluster (stage 1), same as cluster_manta.nf / cluster_wham.nf ---

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

    // --- Cross-caller merge (stage 2) ---

    caller_vcfs = VCFS_CLUSTER_SVCLUSTER_MANTA.out.clustered_vcf
        .mix(VCFS_CLUSTER_SVCLUSTER_WHAM.out.clustered_vcf)

    clustering_config_part1 = channel.of([ [id: 'cc1'], file(params.clustering_config_part1) ])
    stratification_config_part1 = channel.of([ [id: 'sc1'], file(params.stratification_config_part1) ])
    clustering_config_part2 = channel.of([ [id: 'cc2'], file(params.clustering_config_part2) ])
    stratification_config_part2 = channel.of([ [id: 'sc2'], file(params.stratification_config_part2) ])
    track_bed_files = channel.of([
        [id: 'tracks'],
        [
            file(params.clustering_track_sr),
            file(params.clustering_track_sd),
            file(params.clustering_track_rm)
        ]
    ])

    VCFS_COMBINE_BATCHES(
        caller_vcfs,
        ped,
        contig_list,
        fasta,
        fasta_fai,
        dict,
        clustering_config_part1,
        stratification_config_part1,
        clustering_config_part2,
        stratification_config_part2,
        track_bed_files,
        ['SR', 'SD', 'RM'],
        'cohort',
        "${pipelineRoot}/bin/ploidy_table_from_ped.py",
        "${pipelineRoot}/bin/format_gatk_vcf_for_svtk.py"
    )

    VCFS_COMBINE_BATCHES.out.combined_vcf.view { meta, vcf, tbi -> "Combined cross-caller VCF: ${vcf}" }
}
