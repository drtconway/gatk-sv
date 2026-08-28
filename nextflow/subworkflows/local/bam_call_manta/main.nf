//
// Manta germline SV calling, one BAM/CRAM per sample.
//
// GATK-SV's own Manta task (wdl/Manta.wdl) runs one Manta invocation per
// sample and keeps only the diploid SV VCF -- candidate VCFs are Manta
// intermediates GATK-SV never consumes downstream, so we drop them here
// too rather than plumbing them through for no consumer.
//
// GATK-SV's task also rewrites the output VCF's sample column to
// ped_id (wdl/Manta.wdl:152, `bcftools reheader`) since Manta has no
// --sample-name flag and otherwise writes the BAM's own SM tag. This *is*
// implemented here (MANTA_FIX_SAMPLE_ID) -- found via a real KeyError
// downstream (SVCluster's ploidy-table lookup is keyed by the pedigree
// file's individual_id) when it wasn't, on a run against real data. See
// modules/local/manta/fix_sample_id/main.nf for the known gap this
// doesn't yet cover (Manta's convertInversion.py post-processing).
//
// GATK-SV also restricts Manta's own calling to primary contigs (+ mito)
// via --callRegions (wdl/Manta.wdl:137, the manta_region_bed resource).
// Implemented here (target_bed on MANTA_GERMLINE). On top of that, output
// is filtered a second time with VCF_PRIMARY_CONTIGS_ONLY -- callRegions
// constrains Manta's own calling, but doesn't guarantee zero ALT-contig
// records in practice; a real run on real data produced a
// chr1_KI270706v1_random record even with GATK-SV's own equivalent
// upstream restriction. The ploidy table downstream is keyed only by
// primary_contigs_list, so any other contig fails there with a KeyError.
//

include { MANTA_GERMLINE          } from '../../../modules/nf-core/manta/germline/main'
include { MANTA_FIX_SAMPLE_ID     } from '../../../modules/local/manta/fix_sample_id/main'
include { VCF_PRIMARY_CONTIGS_ONLY } from '../../../modules/local/vcf_primary_contigs_only/main'

workflow BAM_CALL_MANTA {
    take:
    bam                // channel: [mandatory] [ meta, bam, bai ]
    fasta              // channel: [mandatory] [ meta2, fasta ]
    fasta_fai          // channel: [mandatory] [ meta3, fai ]
    manta_region_bed     // channel: [mandatory] [ meta4, bed, bed_tbi ] -- primary contigs + mito, for Manta's --callRegions
    primary_contigs_list  // channel: [mandatory] [ meta5, contig_list ] -- for post-hoc filtering

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

    MANTA_FIX_SAMPLE_ID(MANTA_GERMLINE.out.diploid_sv_vcf)

    VCF_PRIMARY_CONTIGS_ONLY(
        MANTA_FIX_SAMPLE_ID.out.vcf,
        primary_contigs_list.map { meta, c -> c }.collect()
    )

    vcf = VCF_PRIMARY_CONTIGS_ONLY.out.vcf
    vcf_tbi = VCF_PRIMARY_CONTIGS_ONLY.out.tbi

    emit:
    vcf       // channel: [ meta, vcf ]
    vcf_tbi   // channel: [ meta, vcf_tbi ]
}
