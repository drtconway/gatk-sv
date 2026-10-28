//
// Reduce a PED file to what the ploidy table needs: individual ID and sex,
// with family, parents and phenotype zeroed and rows sorted by ID. Runs
// before every PLOIDY_TABLE_FROM_PED, which uses `cache 'deep'`
// (conf/modules.config), so edits that don't change anyone's sex
// (parentage, phenotype, family IDs, row order) leave the ploidy tables,
// and everything downstream of them, cached. Without this, any edit to the
// PED reran every cohort-level step (README: "PED edits no longer rerun
// the cohort").
//
process PED_SEXES {
    tag "${meta.id}"
    label 'process_single'

    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/0b/0b4d52ca9a56d07be3f78a12af654e5116f5112908dba277e6796fd9dfb83fe5/data' :
        'community.wave.seqera.io/library/bcftools_htslib:1.23.1--9f08ec665533d64a' }"

    input:
    tuple val(meta), path(ped)

    output:
    tuple val(meta), path("${task.ext.prefix ?: meta.id}.sexes.ped"), emit: ped

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    awk -F '\\t' -v OFS='\\t' '!/^#/ && NF >= 5 { print "0", \$2, "0", "0", \$5, "0" }' ${ped} \\
        | LC_ALL=C sort -t \$'\\t' -k2,2 > ${prefix}.sexes.ped
    """

    // The real command, not a placeholder: -stub-run then exercises the
    // deep-cache behaviour too.
    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    awk -F '\\t' -v OFS='\\t' '!/^#/ && NF >= 5 { print "0", \$2, "0", "0", \$5, "0" }' ${ped} \\
        | LC_ALL=C sort -t \$'\\t' -k2,2 > ${prefix}.sexes.ped
    """
}
