//
// Resolve complex SVs: link inversion, translocation and CNV breakpoints
// that belong to one event into CPX/INV/CTX records. Mirrors GATK-SV's
// ResolveComplexVariants (wdl/ResolveComplexVariants.wdl, with
// wdl/ResolveCpxSv.wdl), the second step of MakeCohortVcf, on
// CombineBatches' output.
//
// Per contig, as upstream:
//   1. RESOLVE_PREP_CONTIG: the contig's inversions only, and the contig's
//      full VCF after the same-breakpoint overlap filter
//   2. Two resolution passes, inversions only (wide preclustering, 50Mb /
//      0.1) and all variants (2kb / ~0), each:
//        SHARD_VCF_CPX       candidate clusters, sharded by variant ID
//        SVTK_RESOLVE_SHARD  svtk resolve per shard
//        CONCAT_RESOLVED_SHARDS
//      A pass with no records (e.g. no inversions) passes its input on.
//   3. INTEGRATE_RESOLVED_VCFS: merge the two passes and rename
// Then genome-wide: concat the contigs (upstream's ReshardVcf, which
// re-splits them by contig for the per-contig steps after) and map the SR
// flag lists onto the new IDs (UPDATE_SR_LISTS, upstream's UpdateSrList).
//
// Every task runs in GATK-SV's sv-pipeline image with its bundled scripts.
// Upstream's per-shard tasks are folded into SVTK_RESOLVE_SHARD to keep
// the Slurm job count down; the commands are the same.
//

include { GET_SE_CUTOFF           } from '../../../modules/local/resolve_complex/get_se_cutoff/main'
include { RESOLVE_PREP_CONTIG     } from '../../../modules/local/resolve_complex/prep_contig/main'
include { SHARD_VCF_CPX as SHARD_VCF_CPX_INV } from '../../../modules/local/resolve_complex/shard_vcf_cpx/main'
include { SHARD_VCF_CPX as SHARD_VCF_CPX_ALL } from '../../../modules/local/resolve_complex/shard_vcf_cpx/main'
include { SVTK_RESOLVE_SHARD      } from '../../../modules/local/resolve_complex/svtk_resolve_shard/main'
include { CONCAT_VCFS as CONCAT_RESOLVED_SHARDS } from '../../../modules/local/concat_vcfs/main'
include { INTEGRATE_RESOLVED_VCFS } from '../../../modules/local/resolve_complex/integrate_resolved_vcfs/main'
include { CONCAT_VCFS as CONCAT_COMPLEX_RESOLVE; contigOrdered } from '../../../modules/local/concat_vcfs/main'
include { UPDATE_SR_LISTS         } from '../../../modules/local/resolve_complex/update_sr_lists/main'

workflow VCF_RESOLVE_COMPLEX_VARIANTS {
    take:
    vcf                  // channel: [ meta, vcf, vcf_tbi ] -- CombineBatches' combined_vcf (svtk format)
    bothside_pass        // channel: [ meta, txt ] -- CombineBatches' BOTHSIDES_SUPPORT ID list
    background_fail      // channel: [ meta, txt ] -- CombineBatches' HIGH_SR_BACKGROUND ID list
    disc_file            // channel: [ meta, pe_txt_gz, tbi ] -- panel-wide PE evidence (vcfs_merge_evidence's merged_pe)
    rf_cutoffs           // channel: [ meta, tsv ] -- ADJUDICATE_SV's cutoffs
    contig_list          // channel: [ meta, contig_list ]
    reference_dict       // channel: [ meta, dict ]
    cytobands            // val: [ bed_gz, tbi ]
    mei_bed              // val: bed_gz
    pe_exclude_list      // val: [ bed_gz, tbi ] -- params.pesr_exclude_list, as upstream
    cohort_name          // val: String, the variant ID prefix
    max_shard_size       // val: Integer, records per svtk resolve shard (upstream: 500)

    main:
    contigs = contig_list.flatMap { meta, c -> file(c).readLines().findAll { it } }

    bothside_value   = bothside_pass.map { meta, f -> f }.collect()
    background_value = background_fail.map { meta, f -> f }.collect()
    disc_value       = disc_file.map { meta, f, tbi -> [ [ f, tbi ] ] }.collect().map { it[0] }
    dict_value       = reference_dict.map { meta, d -> d }.collect()

    GET_SE_CUTOFF(rf_cutoffs)
    se_cutoff = GET_SE_CUTOFF.out.cutoff.map { meta, f -> f.text.trim() }.first()

    // --- 1. per-contig inputs ---
    prep_input = vcf
        .combine(contigs)
        .map { meta, v, tbi, contig -> [ [id: "${cohort_name}.${contig}", contig: contig, cohort: cohort_name], v, tbi ] }

    RESOLVE_PREP_CONTIG(prep_input, bothside_value, background_value)

    // --- 2. the two resolution passes ---
    inv_input = RESOLVE_PREP_CONTIG.out.inv_vcf
        .map { meta, v, tbi -> [ meta + [id: "${meta.id}.inv_only", pass: 'inv'], v, tbi ] }
    all_input = RESOLVE_PREP_CONTIG.out.all_vcf
        .map { meta, v, tbi -> [ meta + [id: "${meta.id}.all", pass: 'all'], v, tbi ] }

    // Preclustering distance/overlap per pass, as upstream.
    SHARD_VCF_CPX_INV(inv_input, 50000000, 0.1, max_shard_size)
    SHARD_VCF_CPX_ALL(all_input, 2000, 0.000000001, max_shard_size)

    pass_inputs = inv_input.mix(all_input).map { meta, v, tbi -> [ meta.id, meta, v, tbi ] }

    // One task per shard; meta.key ties a shard back to its pass. The
    // variant prefix is upstream's ("<cohort>_<pass>_" + "_shard_<i>_");
    // INTEGRATE_RESOLVED_VCFS renames everything afterwards anyway.
    shards = SHARD_VCF_CPX_INV.out.vids.mix(SHARD_VCF_CPX_ALL.out.vids)
        .flatMap { meta, lists ->
            (lists instanceof List ? lists : [ lists ]).collect { l ->
                def i = (l.name =~ /shard_(\d+)\.list$/)[0][1]
                [ meta.id, meta + [id: "${meta.id}.shard_${i}", key: meta.id, variant_prefix: "${meta.cohort}_${meta.pass}__shard_${i}_"], l ]
            }
        }
        .combine(pass_inputs, by: 0)
        .map { key, smeta, vids, pmeta, v, tbi -> [ smeta, v, tbi, vids ] }

    SVTK_RESOLVE_SHARD(shards, disc_value, dict_value, cytobands, mei_bed, pe_exclude_list, se_cutoff)

    // Sorted by name so the concat's inputs (and its hash) are stable.
    concat_input = SVTK_RESOLVE_SHARD.out.vcf
        .map { meta, v, tbi -> [ meta.key, [ v, tbi ] ] }
        .groupTuple()
        .map { key, pairs ->
            def sorted = pairs.sort { it[0].name }
            [ [id: "${key}.resolved", key: key], sorted.collect { it[0] }, sorted.collect { it[1] } ]
        }

    CONCAT_RESOLVED_SHARDS(concat_input)

    // A pass with no shards keeps its input VCF. Indexed, not destructured:
    // an unmatched row's padding isn't guaranteed to be two nulls.
    resolved = pass_inputs
        .join(CONCAT_RESOLVED_SHARDS.out.vcf.map { meta, v, tbi -> [ meta.key, v, tbi ] }, remainder: true)
        .map { row ->
            def (key, pmeta, v, tbi) = row
            row.size() > 4 && row[4] != null ? [ pmeta, row[4], row[5] ] : [ pmeta, v, tbi ]
        }

    // --- 3. integrate the passes and rename, per contig ---
    integrate_input = resolved.filter { meta, v, tbi -> meta.pass == 'inv' }.map { meta, v, tbi -> [ meta.contig, v, tbi ] }
        .join(resolved.filter { meta, v, tbi -> meta.pass == 'all' }.map { meta, v, tbi -> [ meta.contig, v, tbi ] })
        .map { contig, iv, itbi, av, atbi -> [ [id: "${cohort_name}.${contig}.renamed", contig: contig, cohort: cohort_name], iv, itbi, av, atbi ] }

    INTEGRATE_RESOLVED_VCFS(integrate_input)

    // --- genome-wide ---
    CONCAT_COMPLEX_RESOLVE(
        INTEGRATE_RESOLVED_VCFS.out.vcf
            .map { meta, v, tbi -> [ meta.contig, v, tbi ] }
            .collect(flat: false)
            .map { rows -> [ rows ] }
            .combine(contig_list.map { meta, c -> [ file(c).readLines().findAll { it } ] })
            .map { rows, names -> contigOrdered([id: "${cohort_name}.complex_resolve"], rows, names) }
    )

    UPDATE_SR_LISTS(CONCAT_COMPLEX_RESOLVE.out.vcf, bothside_value, background_value)

    emit:
    vcf             = CONCAT_COMPLEX_RESOLVE.out.vcf          // channel: [ meta, vcf, vcf_tbi ]
    bothside_pass   = UPDATE_SR_LISTS.out.bothside_pass       // channel: [ meta, txt ]
    background_fail = UPDATE_SR_LISTS.out.background_fail     // channel: [ meta, txt ]
    dropped_vcfs    = RESOLVE_PREP_CONTIG.out.dropped_vcf     // channel: [ meta, vcf ] -- per contig, records the breakpoint overlap filter dropped
}
