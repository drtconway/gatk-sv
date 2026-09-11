//
// Cross-check the sample sheet (already parsed into per-sample [meta, bam,
// bai] tuples by utils_input_channels) against the pedigree file, and fail
// fast with every mismatch listed at once if they disagree -- rather than
// letting a mismatch surface many stages later as an opaque KeyError deep
// in some ped_id-keyed lookup (ploidy tables, VCF sample-column rewrites,
// --sample-name on evidence-collection walkers -- see CLAUDE.md's own
// ped_id-vs-sample_id note). Several real HPC runs hit exactly this class
// of bug (a sample sheet row's ped_id with no matching PED individual_id,
// or vice versa) with no validation catching it before the pipeline had
// already spent hours on earlier stages.
//
// Deliberately narrow: pure text cross-referencing between the two input
// files, no BAM/CRAM header access (e.g. verifying a BAM's own @RG SM: tag
// matches its assigned sample identity is a different, heavier class of
// check -- real upstream data mislabeling, not sample-sheet/PED
// inconsistency -- and not what this subworkflow does).
//
// Not every workflows/*.nf entry point needs this: only those that
// actually consume the PED file downstream (cluster_manta.nf,
// cluster_wham.nf, combine_batches.nf, generate_batch_metrics.nf,
// genotype_batch.nf) -- so this is a separate subworkflow from
// utils_input_channels (used by all ten entry points, including PED-free
// ones like call_manta.nf), called once, right after it.
//

workflow VALIDATE_SAMPLESHEET_PED {
    take:
    samples   // channel: [ meta, bam, bai ] -- from utils_input_channels, meta has .id (sample_id) and .ped_id
    ped_path  // val: params.ped

    main:
    def ped_file = file(ped_path)

    // individual_id -> [paternal_id, maternal_id], skipping the leading
    // '#'-commented header line PED files conventionally carry (see
    // README.md#pedigree-file's own example) and any other comment/blank
    // lines.
    def ped_individuals = [:]
    ped_file.readLines().each { line ->
        if (!line || line.startsWith('#')) return
        def fields = line.split('\t')
        if (fields.size() < 6) {
            error "validate_samplesheet_ped: malformed PED line (expected 6 tab-separated fields, got ${fields.size()}) in ${ped_path}: '${line}'"
        }
        def (fam, individual_id, paternal_id, maternal_id, sex, phenotype) = fields
        ped_individuals[individual_id] = [paternal_id, maternal_id]
    }

    // Collect every sample sheet row's (sample_id, effective ped_id) pair,
    // then do every cross-check at once outside the channel-processing
    // machinery -- this needs the *complete* sample list, not a per-row
    // check, since "every PED individual_id maps back to exactly one
    // sample sheet row" and similar are whole-list properties.
    samples
        .map { meta, bam, bai -> [ sample_id: meta.id, ped_id: meta.ped_id ] }
        .collect()
        .map { rows ->
            def errors = []

            // 1. Every sample sheet row's ped_id must exist in the PED file.
            rows.each { row ->
                if (!ped_individuals.containsKey(row.ped_id)) {
                    errors << "sample sheet row '${row.sample_id}' has ped_id '${row.ped_id}', which is not a PED individual_id in ${ped_path}"
                }
            }

            // 2. Every PED individual_id referenced by the sample sheet
            // should map back to exactly one row -- catches two different
            // sample_ids accidentally sharing one ped_id (the exact
            // failure mode a real HPC run hit: one sample's VCF ended up
            // with another sample's genotype column because both
            // resolved to the same ped_id downstream).
            def ped_id_to_samples = rows.groupBy { it.ped_id }
            ped_id_to_samples.each { ped_id, group ->
                if (group.size() > 1) {
                    def sample_ids = group.collect { it.sample_id }.join(', ')
                    errors << "ped_id '${ped_id}' is used by more than one sample sheet row: ${sample_ids}"
                }
            }

            // 3. PED parent references that aren't '0' should themselves
            // be a real individual_id in the same file -- a dangling
            // paternal_id/maternal_id doesn't break the sample-sheet link
            // directly, but breaks trio logic downstream just as opaquely.
            ped_individuals.each { individual_id, parents ->
                def (paternal_id, maternal_id) = parents
                if (paternal_id != '0' && !ped_individuals.containsKey(paternal_id)) {
                    errors << "PED individual '${individual_id}' has paternal_id '${paternal_id}', which is not itself a PED individual_id in ${ped_path}"
                }
                if (maternal_id != '0' && !ped_individuals.containsKey(maternal_id)) {
                    errors << "PED individual '${individual_id}' has maternal_id '${maternal_id}', which is not itself a PED individual_id in ${ped_path}"
                }
            }

            if (errors) {
                error "Sample sheet / PED file inconsistency (${ped_path}):\n" + errors.collect { "  - ${it}" }.join('\n')
            }
        }
}
