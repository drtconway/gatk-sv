#!/usr/bin/env nextflow

//
// Standalone entry point: sample sheet -> Wham VCFs (contig-restricted,
// sample ID/TAGS fixed to match the pipeline's sample_id).
//
// Run directly, e.g.:
//   nextflow run workflows/call_wham.nf -profile apptainer \
//     --input samplesheet.tsv --fasta ref.fa --fasta_fai ref.fa.fai \
//     --primary_contigs_list primary_contigs.list
//
// Same composable-caller pattern as call_manta.nf -- see
// ../README.md#composability.
//

include { samplesheetToList } from 'plugin/nf-schema'
include { BAM_CALL_WHAM     } from '../subworkflows/local/bam_call_wham/main'

workflow {
    main:
    // projectDir resolves to the directory of *this* script (workflows/)
    // when run standalone, so the pipeline root -- where assets/ lives --
    // has to be derived rather than assumed to be projectDir itself.
    def pipelineRoot = "${projectDir}/.."

    if (!params.primary_contigs_list) {
        error "call_wham.nf requires --primary_contigs_list (see conf/modules.config)"
    }

    // 'family_id' and 'sample_id' are both tagged meta: in the schema, so
    // samplesheetToList already merges them into one meta map per row:
    // [ [family_id:.., id:..], bam, bai ]
    samples = channel
        .fromList(samplesheetToList(params.input, "${pipelineRoot}/assets/schema_samplesheet.json"))
        .map { meta, bam, bai -> [ meta, file(bam), file(bai) ] }

    fasta     = channel.of([ [ id: 'reference' ], file(params.fasta) ])
    fasta_fai = channel.of([ [ id: 'reference' ], file(params.fasta_fai) ])

    BAM_CALL_WHAM(
        samples,
        fasta,
        fasta_fai
    )

    BAM_CALL_WHAM.out.vcf.view { meta, vcf -> "Wham VCF for ${meta.id}: ${vcf}" }
}
