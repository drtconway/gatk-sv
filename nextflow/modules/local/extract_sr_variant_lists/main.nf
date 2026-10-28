//
// List the variant IDs carrying INFO/HIGH_SR_BACKGROUND and
// INFO/BOTHSIDES_SUPPORT. Reproduces wdl/CombineBatches.wdl's
// ExtractSRVariantLists task; ResolveComplexVariants and CleanVcf take the
// two lists as inputs. A flag missing from the header gives an empty list
// (with a warning) rather than a bcftools error.
//
process EXTRACT_SR_VARIANT_LISTS {
    tag "${meta.id}"
    label 'process_single'

    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/0b/0b4d52ca9a56d07be3f78a12af654e5116f5112908dba277e6796fd9dfb83fe5/data' :
        'community.wave.seqera.io/library/bcftools_htslib:1.23.1--9f08ec665533d64a' }"

    input:
    tuple val(meta), path(vcf), path(vcf_index)

    output:
    tuple val(meta), path("*.high_sr_background.txt"), emit: background_fail
    tuple val(meta), path("*.bothsides_sr_support.txt"), emit: bothside_pass

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    // Header captured first: `bcftools view -h | grep -q` would SIGPIPE
    // bcftools under pipefail.
    """
    hdr=\$(bcftools view -h ${vcf})
    fields='%ID'
    for flag in HIGH_SR_BACKGROUND BOTHSIDES_SUPPORT; do
        if grep -q "^##INFO=<ID=\${flag}," <<< "\$hdr"; then
            fields="\${fields}\\t%INFO/\${flag}"
        else
            echo "WARNING: INFO/\${flag} not in ${vcf}'s header; its list will be empty" >&2
            fields="\${fields}\\t."
        fi
    done
    bcftools query -f "\${fields}\\n" ${vcf} > flags.txt
    awk -F '\\t' '\$2 != "." { print \$1 }' flags.txt > ${prefix}.high_sr_background.txt
    awk -F '\\t' '\$3 != "." { print \$1 }' flags.txt > ${prefix}.bothsides_sr_support.txt
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.high_sr_background.txt ${prefix}.bothsides_sr_support.txt
    """
}
