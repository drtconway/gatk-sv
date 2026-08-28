//
// Manta germline SV calling, one BAM/CRAM per sample.
//
// GATK-SV's own Manta task (wdl/Manta.wdl) runs one Manta invocation per
// sample and keeps only the diploid SV VCF -- candidate VCFs are Manta
// intermediates GATK-SV never consumes downstream, so we drop them here
// too rather than plumbing them through for no consumer.
//
// GATK-SV also restricts Manta's own calling to primary contigs (+ mito)
// via --callRegions (wdl/Manta.wdl:137, the manta_region_bed resource).
// Implemented here (target_bed on MANTA_GERMLINE).
//
// Manta's raw output is then standardized with SVTK_STANDARDIZE (wraps
// `svtk standardize`, GATK-SV's own StandardizeVCFs task,
// wdl/PESRPreprocessing.wdl) -- this single step does the sample-column
// rewrite to ped_id, primary-contig restriction, min-size filtering, and
// sets standard INFO fields (including ALGORITHMS, which SVCluster
// requires downstream). See modules/local/svtk_standardize/main.nf: this
// replaced three separate ad-hoc reconstructions (bcftools-reheader
// rename, a post-hoc contig filter) that were each independently found
// incomplete via real KeyErrors/GATKExceptions on real data -- none of
// them set ALGORITHMS, which SVCluster requires.
//
// Not implemented: GATK-SV also runs Manta's raw output through
// convertInversion.py (converts Manta's inversion-signaling BND pairs
// into proper INV records) before standardization. See README known gaps.
//

include { MANTA_GERMLINE  } from '../../../modules/nf-core/manta/germline/main'
include { SVTK_STANDARDIZE } from '../../../modules/local/svtk_standardize/main'

workflow BAM_CALL_MANTA {
    take:
    bam                // channel: [mandatory] [ meta, bam, bai ]
    fasta              // channel: [mandatory] [ meta2, fasta ]
    fasta_fai          // channel: [mandatory] [ meta3, fai ]
    manta_region_bed   // channel: [mandatory] [ meta4, bed, bed_tbi ] -- primary contigs + mito, for Manta's --callRegions
    primary_contigs_fai // channel: [mandatory] [ meta5, contigs_fai ] -- for svtk standardize --contigs
    min_size             // val: Int, minimum SV size to retain (GATK-SV default: 50)

    main:
    manta_input = bam
        .combine(manta_region_bed.map { meta, bed, tbi -> [ bed, tbi ] })
        .map { meta, bam_file, bai_file, bed, tbi ->
            [ meta, bam_file, bai_file, bed, tbi ]
        }

    // fasta/fasta_fai are one reference shared across every sample, but
    // MANTA_GERMLINE declares them as ordinary (non-value) channel inputs.
    // Without .collect() they're each consumed once and then exhausted, so
    // Nextflow silently pairs positionally and drops every sample after
    // the first -- .collect() turns them into a single broadcastable value
    // reused for every item in manta_input.
    MANTA_GERMLINE(
        manta_input,
        fasta.collect(),
        fasta_fai.collect(),
        []  // no config override
    )
    // NB: MANTA_GERMLINE.out.versions_manta is deliberately not consumed.
    // It's a `topic: versions, emit: ..., eval("configManta.py --version")`
    // output -- Nextflow's `eval` output type always shells out for real,
    // even under -stub-run, so referencing it breaks stub-mode validation
    // when the real tool isn't installed. Wire it up once we're running
    // against real containers and actually need version tracking.

    SVTK_STANDARDIZE(
        MANTA_GERMLINE.out.diploid_sv_vcf,
        primary_contigs_fai.map { meta, c -> c }.collect(),
        'manta',
        min_size
    )

    vcf = SVTK_STANDARDIZE.out.vcf
    vcf_tbi = SVTK_STANDARDIZE.out.tbi

    emit:
    vcf       // channel: [ meta, vcf ]
    vcf_tbi   // channel: [ meta, vcf_tbi ]
}
