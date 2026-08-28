#!/usr/bin/env nextflow

//
// Standalone entry point: sample sheet -> Manta VCFs -> one clustered
// site VCF across all samples (stage 1 harmonisation -- see
// subworkflows/local/vcfs_cluster_svcluster).
//
// Run directly, e.g.:
//   nextflow run workflows/cluster_manta.nf -profile apptainer \
//     --input samplesheet.tsv --fasta ref.fa --fasta_fai ref.fa.fai \
//     --primary_contigs_list primary_contigs.list \
//     --primary_contigs_fai contigs.fai --ped cohort.ped \
//     --pesr_exclude_intervals exclude.bed.gz \
//     --pesr_exclude_intervals_tbi exclude.bed.gz.tbi \
//     --manta_region_bed manta_region.bed.gz \
//     --manta_region_bed_tbi manta_region.bed.gz.tbi
//

include { UTILS_INPUT_CHANNELS  } from '../subworkflows/local/utils_input_channels/main'
include { BAM_CALL_MANTA        } from '../subworkflows/local/bam_call_manta/main'
include { VCFS_CLUSTER_SVCLUSTER } from '../subworkflows/local/vcfs_cluster_svcluster/main'

workflow {
    main:
    def pipelineRoot = "${projectDir}/.."

    [
        'primary_contigs_list', 'primary_contigs_fai', 'ped',
        'pesr_exclude_intervals', 'pesr_exclude_intervals_tbi',
        'reference_dict', 'manta_region_bed', 'manta_region_bed_tbi'
    ].each { p ->
        if (!params[p]) {
            error "cluster_manta.nf requires --${p}"
        }
    }

    UTILS_INPUT_CHANNELS(params.input, pipelineRoot, params.fasta, params.fasta_fai)

    // primary_contigs_list (plain list, for vcfs_cluster_svcluster's
    // ploidy-table/format-conversion scripts) and primary_contigs_fai
    // (.fai format, for svtk standardize) are both needed here.
    contig_list = channel.of([ [id: 'contigs'], file(params.primary_contigs_list) ])
    contigs_fai = channel.of([ [id: 'contigs'], file(params.primary_contigs_fai) ])
    manta_region_bed = channel.of([
        [id: 'manta_region'],
        file(params.manta_region_bed),
        file(params.manta_region_bed_tbi)
    ])

    BAM_CALL_MANTA(
        UTILS_INPUT_CHANNELS.out.samples,
        UTILS_INPUT_CHANNELS.out.fasta,
        UTILS_INPUT_CHANNELS.out.fasta_fai,
        manta_region_bed,
        contigs_fai,
        params.min_svsize
    )

    ped = channel.of([ [id: 'ped'], file(params.ped) ])
    exclude_intervals = channel.of([
        [id: 'exclude'],
        file(params.pesr_exclude_intervals),
        file(params.pesr_exclude_intervals_tbi)
    ])
    dict = channel.of([ [id: 'reference'], file(params.reference_dict) ])

    VCFS_CLUSTER_SVCLUSTER(
        BAM_CALL_MANTA.out.vcf,
        ped,
        contig_list,
        exclude_intervals,
        UTILS_INPUT_CHANNELS.out.fasta,
        UTILS_INPUT_CHANNELS.out.fasta_fai,
        dict,
        'manta',
        params.min_svsize,
        "${pipelineRoot}/bin/ploidy_table_from_ped.py",
        "${pipelineRoot}/bin/format_svtk_vcf_for_gatk.py",
        "${pipelineRoot}/bin/format_gatk_vcf_for_svtk.py"
    )

    VCFS_CLUSTER_SVCLUSTER.out.clustered_vcf.view { meta, vcf, tbi -> "Clustered Manta VCF: ${vcf}" }
}
