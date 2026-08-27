//
// Wham germline SV calling, one BAM/CRAM per sample.
//
// GATK-SV's own Wham task (wdl/Whamg.wdl) does more than a plain whamg
// invocation:
//   1. Restricts calling to primary_contigs_list via `-c` (for read-stat
//      estimation, not correctness-critical region exclusion). Implemented
//      here via task.ext.args in conf/modules.config, which reads
//      params.primary_contigs_list -- see the note there for why this
//      couldn't be wired through as an ordinary process input.
//   2. Scatters over an include_bed_file whitelist of regions, then
//      concatenates -- this bounds runtime and avoids regions Wham
//      struggles with. NOT implemented here yet; we run whamg genome-wide
//      in one shot. Revisit if runtime on real data warrants it.
//   3. Rewrites the output VCF's sample column and TAGS INFO field to the
//      pipeline's sample_id, since Wham defaults to the BAM's own SM tag
//      for both. This *is* implemented here (see WHAM_FIX_SAMPLE_ID) --
//      GATK-SV's downstream svtk standardize_vcf reads TAGS specifically
//      for WHAM VCFs, so a mismatch here silently misattributes calls.
//

include { WHAMG              } from '../../../modules/nf-core/whamg/main'
include { WHAM_FIX_SAMPLE_ID } from '../../../modules/local/wham/fix_sample_id/main'

workflow BAM_CALL_WHAM {
    take:
    bam        // channel: [mandatory] [ meta, bam, bai ]
    fasta      // channel: [mandatory] [ meta2, fasta ]
    fasta_fai  // channel: [mandatory] [ meta3, fai ]

    main:
    // fasta/fasta_fai are one reference shared across every sample; see
    // the same .collect() note in bam_call_manta/main.nf for why this is
    // needed (ordinary channel inputs get exhausted after one sample
    // otherwise). WHAMG also wants bare paths, not [meta, path] tuples, so
    // strip the meta *before* collecting -- collect() first would wrap the
    // whole [meta, fa] tuple in a one-element list instead.
    WHAMG(
        bam,
        fasta.map { meta, fa -> fa }.collect(),
        fasta_fai.map { meta, fai -> fai }.collect()
    )

    fix_input = WHAMG.out.vcf.join(WHAMG.out.tbi)

    WHAM_FIX_SAMPLE_ID(fix_input)

    vcf = WHAM_FIX_SAMPLE_ID.out.vcf
    vcf_tbi = WHAM_FIX_SAMPLE_ID.out.tbi

    emit:
    vcf      // channel: [ meta, vcf ]
    vcf_tbi  // channel: [ meta, vcf_tbi ]
}
