//
// Map the SR flag lists' variant IDs onto the resolved VCF's new IDs (via
// INFO/MEMBERS). Reproduces TasksMakeCohortVcf.wdl's UpdateSrList, run for
// both lists. Run once on the genome-wide VCF, not per contig as upstream:
// IDs not in the VCF are kept as they are, so per-contig runs over these
// genome-wide lists would repeat every ID once per contig.
//
process UPDATE_SR_LISTS {
    tag "${meta.id}"
    label 'process_single'

    container 'us.gcr.io/broad-dsde-methods/gatk-sv/sv-pipeline:2026-06-02-v1.1-24f0fd08'

    input:
    tuple val(meta), path(vcf), path(vcf_index)
    path(bothside_pass_list)
    path(background_fail_list)

    output:
    tuple val(meta), path("${task.ext.prefix ?: meta.id}.sr_bothside_pass.updated3.txt"), emit: bothside_pass
    tuple val(meta), path("${task.ext.prefix ?: meta.id}.sr_background_fail.updated3.txt"), emit: background_fail

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    // gawk (ARGIND), as in the sv-pipeline image.
    """
    svtk vcf2bed ${vcf} int.bed -i MEMBERS --no-samples --no-header

    update() {
        awk -F'[,\\t]' -v OFS='\\t' \\
            '{
                if (ARGIND==1) for (i=6; i<=NF; ++i) MAP[\$i]=\$4;
                else if (\$NF in MAP) print \$0,MAP[\$NF];
                else print \$0,\$NF;
            }' int.bed "\$1" \\
            | sort -T . -k1,1n > "\$2"
    }
    update ${bothside_pass_list} ${prefix}.sr_bothside_pass.updated3.txt
    update ${background_fail_list} ${prefix}.sr_background_fail.updated3.txt
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.sr_bothside_pass.updated3.txt ${prefix}.sr_background_fail.updated3.txt
    """
}
