#!/usr/bin/env nextflow

//
// Standalone entry point: sample sheet -> Manta diploid SV VCFs.
//
// Run directly, e.g.:
//   nextflow run workflows/call_manta.nf -profile apptainer \
//     --input samplesheet.tsv --fasta ref.fa --fasta_fai ref.fa.fai
//
// This is the first slice of the composable pipeline described in
// ../README.md#composability: no panel bundle, harmonisation, or
// genotyping yet -- just the caller subworkflow wired to a sample sheet.
//

include { samplesheetToList } from 'plugin/nf-schema'
include { BAM_CALL_MANTA    } from '../subworkflows/local/bam_call_manta/main'

workflow {
    main:
    // projectDir resolves to the directory of *this* script (workflows/)
    // when run standalone, so the pipeline root -- where assets/ lives --
    // has to be derived rather than assumed to be projectDir itself.
    def pipelineRoot = "${projectDir}/.."

    // 'family_id' and 'sample_id' are both tagged meta: in the schema, so
    // samplesheetToList already merges them into one meta map per row:
    // [ [family_id:.., id:..], bam, bai ]
    samples = channel
        .fromList(samplesheetToList(params.input, "${pipelineRoot}/assets/schema_samplesheet.json"))
        .map { meta, bam, bai -> [ meta, file(bam), file(bai) ] }

    fasta     = channel.of([ [ id: 'reference' ], file(params.fasta) ])
    fasta_fai = channel.of([ [ id: 'reference' ], file(params.fasta_fai) ])

    BAM_CALL_MANTA(
        samples,
        fasta,
        fasta_fai
    )

    BAM_CALL_MANTA.out.vcf.view { meta, vcf -> "Manta diploid SV VCF for ${meta.id}: ${vcf}" }
}
