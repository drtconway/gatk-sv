//
// Wham germline SV calling, one BAM/CRAM per sample.
//
// GATK-SV's own Wham task (wdl/Whamg.wdl) restricts calling to
// primary_contigs_list via `-c` (for read-stat estimation, not
// correctness-critical region exclusion). Implemented here via
// task.ext.args in conf/modules.config, which reads
// params.primary_contigs_list -- see the note there for why this
// couldn't be wired through as an ordinary process input.
//
// NOT implemented: GATK-SV also scatters over an include_bed_file
// whitelist of regions, then concatenates -- this bounds runtime and
// avoids regions Wham struggles with. We run whamg genome-wide in one
// shot instead. Revisit if runtime on real data warrants it.
//
// Wham's raw output is fed directly (unrenamed) into SVTK_STANDARDIZE
// (wraps `svtk standardize`, see modules/local/svtk_standardize/main.nf),
// which does the sample-column rewrite, primary-contig restriction,
// min-size filtering, and sets ALGORITHMS -- same as Manta, see
// bam_call_manta/main.nf.
//
// Wham needs its *raw*, un-renamed output here specifically: Wham writes
// the BAM's own SM tag as both the VCF sample column and the TAGS INFO
// field, and svtk's own Wham standardizer (src/svtk/svtk/standardize/
// std_wham.py) determines genotypes by checking whether each *raw*
// sample name appears in TAGS, then maps to the --sample-names value
// (ped_id) as it writes output -- it does this matching itself; the
// input's sample name and TAGS value need to still agree with each
// other, not with ped_id. An earlier version of this subworkflow
// rewrote the sample column and TAGS to ped_id *before* svtk
// standardize (modules/local/wham/fix_sample_id, now removed) --
// verified unnecessary and the wrong step to do this at, once
// SVTK_STANDARDIZE replaced it: svtk standardize's own TAGS handling
// already does the right thing given Wham's true raw output.
//

include { WHAMG            } from '../../../modules/nf-core/whamg/main'
include { SVTK_STANDARDIZE } from '../../../modules/local/svtk_standardize/main'

workflow BAM_CALL_WHAM {
    take:
    bam                  // channel: [mandatory] [ meta, bam, bai ]
    fasta                // channel: [mandatory] [ meta2, fasta ]
    fasta_fai            // channel: [mandatory] [ meta3, fai ]
    primary_contigs_fai  // channel: [mandatory] [ meta4, contigs_fai ] -- for svtk standardize --contigs
    min_size               // val: Int, minimum SV size to retain (GATK-SV default: 50)

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

    SVTK_STANDARDIZE(
        WHAMG.out.vcf,
        primary_contigs_fai.map { meta, c -> c }.collect(),
        'wham',
        min_size
    )

    vcf = SVTK_STANDARDIZE.out.vcf
    vcf_tbi = SVTK_STANDARDIZE.out.tbi

    emit:
    vcf      // channel: [ meta, vcf ]
    vcf_tbi  // channel: [ meta, vcf_tbi ]
}
