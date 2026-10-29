//
// Resolve one shard of candidate complex SVs. Reproduces, in one task,
// wdl/ResolveCpxSv.wdl's per-shard chain:
//   PullVcfShard          the shard's records, by variant ID
//   ResolvePrep           PE evidence within 2kb of the shard's inversion
//                         breakpoints (PrintSVEvidence), plus the dummy
//                         pairs svtk needs to not break on an empty file
//   SvtkResolve           svtk resolve
//   RestoreUnresolvedCnv  unresolved DEL/DUP back in as plain CNVs, other
//                         unresolved records kept (inversion single-enders
//                         as BND), then postCPX_cleanup.py
// Sorts spill to the work dir, not the node's /tmp (README: "Then
// FILTER_ANNOTATE_VCF (manta) filled the node's /tmp").
//
process SVTK_RESOLVE_SHARD {
    tag "${meta.id}"
    label 'process_single'

    container 'us.gcr.io/broad-dsde-methods/gatk-sv/sv-pipeline:2026-06-02-v1.1-24f0fd08'

    input:
    tuple val(meta), path(vcf), path(vcf_index), path(vids)  // meta.contig, meta.variant_prefix
    tuple path(disc_file), path(disc_file_index)
    path(reference_dict)
    tuple path(cytobands), path(cytobands_index)
    path(mei_bed)
    tuple path(pe_exclude_list), path(pe_exclude_list_index)
    val(se_cutoff)

    output:
    tuple val(meta), path("${task.ext.prefix ?: meta.id}.vcf.gz"), path("${task.ext.prefix ?: meta.id}.vcf.gz.tbi"), emit: vcf

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def chrom = meta.contig
    def avail_mem = task.memory ? (task.memory.mega * 0.8).intValue() : 3072
    """
    mkdir -p tmp

    # --- PullVcfShard ---
    bcftools view --no-version --include ID=@${vids} ${vcf} -Oz -o shard.vcf.gz
    tabix shard.vcf.gz

    # --- ResolvePrep ---
    svtk vcf2bed shard.vcf.gz input.bed --no-samples --no-header
    (grep INV input.bed || true) \\
        | awk -v OFS='\\t' -v buffer=2000 '{ print \$1, \$2-buffer, \$2+buffer"\\n"\$1, \$3-buffer, \$3+buffer }' \\
        | awk -v OFS='\\t' '{ if (\$2<1) \$2=1; print \$1, \$2, \$3 }' \\
        | sort -T ./tmp -Vk1,1 -k2,2n -k3,3n \\
        | bedtools merge -i - \\
        > regions.bed

    if [ -s regions.bed ]; then
        java -Xmx${avail_mem}M -jar \$GATK_JAR PrintSVEvidence \\
            --sequence-dictionary ${reference_dict} \\
            --evidence-file ${disc_file} \\
            -L regions.bed \\
            -O disc.PE.txt
        awk '{ if (\$1==\$4 && \$3==\$6) print }' disc.PE.txt > disc.filtered.PE.txt
    else
        : > disc.filtered.PE.txt
    fi
    # Upstream adds both dummies, including the odd "chr<chrom>" one.
    {
        cat disc.filtered.PE.txt
        printf '%s\\t1\\t+\\t%s\\t2\\t+\\tDUMMY_SAMPLE_IGNORE\\n' ${chrom} ${chrom}
        printf 'chr%s\\t1\\t+\\t%s\\t2\\t+\\tDUMMY_SAMPLE_IGNORE\\n' ${chrom} ${chrom}
    } \\
        | sort -T ./tmp -Vk1,1 -k2,2n -k5,5n -k7,7 \\
        | bgzip -c \\
        > discfile.PE.txt.gz
    tabix -0 -s 1 -b 2 -e 2 -f discfile.PE.txt.gz

    # --- SvtkResolve ---
    svtk resolve \\
        shard.vcf.gz \\
        resolved.vcf \\
        -p ${meta.variant_prefix} \\
        -u unresolved.vcf \\
        -t ./tmp/ \\
        --discfile discfile.PE.txt.gz \\
        --mei-bed ${mei_bed} \\
        --cytobands ${cytobands} \\
        --min-rescan-pe-support ${se_cutoff} \\
        -x ${pe_exclude_list}

    # --- RestoreUnresolvedCnv ---
    bgzip resolved.vcf
    bcftools view --no-version --no-header unresolved.vcf -Oz -o unresolved_records.vcf.gz
    mv resolved.vcf.gz resolved_plus_cnv.tmp.gz

    # Unresolved CNVs, unresolved status wiped
    zcat unresolved_records.vcf.gz \\
        | (grep -F -e "<DEL>" -e "<DUP>" -e "SVTYPE=DEL" -e "SVTYPE=DUP" -e "SVTYPE=CNV" -e "SVTYPE=MCNV" || true) \\
        | sed -r -e 's/;EVENT=[^;]*;/;/' -e 's/;UNRESOLVED[^;]*;/;/g' \\
        | sed -r -e 's/;UNRESOLVED_TYPE[^;]*;/;/g' -e 's/;UNRESOLVED_TYPE[^\\t]*\\t/\\t/g' \\
        | bgzip \\
        >> resolved_plus_cnv.tmp.gz
    # Other unresolved records, unresolved status kept (except inversion single-enders)
    zcat unresolved_records.vcf.gz \\
        | (grep -F -v -e "<DEL>" -e "<DUP>" -e "SVTYPE=DEL" -e "SVTYPE=DUP" -e "SVTYPE=CNV" -e "SVTYPE=MCNV" \\
                     -e "INVERSION_SINGLE_ENDER" || true) \\
        | bgzip \\
        >> resolved_plus_cnv.tmp.gz
    # Inversion single-enders, as BND
    zcat unresolved_records.vcf.gz \\
        | (grep -F -v -e "<DEL>" -e "<DUP>" -e "SVTYPE=DEL" -e "SVTYPE=DUP" -e "SVTYPE=CNV" -e "SVTYPE=MCNV" || true) \\
        | (grep -F -e "INVERSION_SINGLE_ENDER" || true) \\
        | sed -e 's/SVTYPE=INV/SVTYPE=BND/g' \\
        | sed -e 's/END=\\([0-9]*\\)/END=\\1;END2=\\1/' \\
        | bgzip \\
        >> resolved_plus_cnv.tmp.gz

    zcat resolved_plus_cnv.tmp.gz \\
        | vcf-sort -c -t ./tmp \\
        | python /opt/sv-pipeline/04_variant_resolution/scripts/postCPX_cleanup.py /dev/stdin /dev/stdout \\
        | bgzip \\
        > ${prefix}.vcf.gz
    tabix ${prefix}.vcf.gz
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo "" | gzip > ${prefix}.vcf.gz
    touch ${prefix}.vcf.gz.tbi
    """
}
