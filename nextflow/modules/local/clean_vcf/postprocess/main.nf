//
// The rest of CleanVcf for one contig, in upstream's order
// (wdl/CleanVcfChromosome.wdl), unsharded:
//   CleanVcfPostprocess           cleanvcf_postprocess.py, drop working INFO fields
//   DropRedundantCnvs + SortVcf   resolve_cpx_cnv_redundancies.py
//   StitchFragmentedCnvs          StitchFragmentedCNVs.jar, two passes
//   RescueMobileElementDeletions  short +- BNDs over LINE1/HERVK -> DEL:ME
//   AddHighFDRFilters             HIGH_ALGORITHM_FDR on Wham-only DELs
//   AddRetroDelFilters            add_retro_del_filters.py (intron deletions)
//   FinalCleanup                  rename, drop CNVs < 5kb and INFO/MEMBERS
//   FormatVcfToOutput             back to GATK format
// The inline Python is upstream's, with paths passed as arguments.
//
process CLEANVCF_POSTPROCESS {
    tag "${meta.id}"
    label 'process_single'

    container 'us.gcr.io/broad-dsde-methods/gatk-sv/sv-pipeline:2026-06-02-v1.1-24f0fd08'

    input:
    tuple val(meta), path(vcf), path(vcf_index)  // meta.contig, meta.cohort
    path(ploidy_table)
    path(bothside_pass_list)
    path(background_fail_list)
    path(ped_file)          // PED_SEXES' ID+sex table: the script reads only sex
    tuple path(line1_reference), path(hervk_reference), path(intron_reference)
    path(retro_del_script)   // bin/add_retro_del_filters.py (newer than the image's, see bin/README.md)

    output:
    tuple val(meta), path("${task.ext.prefix ?: meta.id}.vcf.gz"), path("${task.ext.prefix ?: meta.id}.vcf.gz.tbi"), emit: vcf

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    // FinalCleanup's ID prefix, as upstream (its task prefix), giving IDs
    // like cohort.chr1.final_cleanup_DEL_chr1_1
    def id_prefix = "${meta.cohort}.${meta.contig}.final_cleanup"
    def avail_mem = task.memory ? (task.memory.mega * 0.8).intValue() : 3072
    """
    mkdir -p tmp

    # --- CleanVcfPostprocess ---
    python /opt/sv-pipeline/04_variant_resolution/scripts/cleanvcf_postprocess.py \\
        -V ${vcf} \\
        -O processed.vcf.gz \\
        --ped-file ${ped_file}
    bcftools annotate -x INFO/MULTIALLELIC,INFO/UNRESOLVED,INFO/EVENT,INFO/REVISED_EVENT,INFO/MULTI_CNV,INFO/varGQ \\
        processed.vcf.gz -o processed.annotated.vcf.gz -O z
    bcftools view -h processed.annotated.vcf.gz | grep "^##" \\
        | grep -v -E "CIPOS|CIEND|RMSSTD|source|bcftools|GATKCommandLine|##FORMAT=<ID=EV>|##ALT=<ID=UNR>|##INFO=<ID=(MULTIALLELIC|UNRESOLVED|EVENT|REVISED_EVENT|MULTI_CNV|varGQ)" \\
        > temp_header.txt
    echo '##INFO=<ID=UNRESOLVED_TYPE,Number=1,Type=String,Description="Class of unresolved variant.">' >> temp_header.txt
    echo '##ALT=<ID=CNV,Description="Copy Number Polymorphism">' >> temp_header.txt
    bcftools view -h processed.annotated.vcf.gz | grep "^#CHROM" > chrom_header.txt
    cat temp_header.txt chrom_header.txt > header.txt
    bcftools reheader -h header.txt processed.annotated.vcf.gz -o postprocess.vcf.gz
    tabix -p vcf postprocess.vcf.gz

    # --- DropRedundantCnvs + SortVcf ---
    python /opt/sv-pipeline/04_variant_resolution/scripts/resolve_cpx_cnv_redundancies.py \\
        postprocess.vcf.gz drop_redundant_cnvs.vcf.gz --temp-dir ./tmp
    bcftools sort --temp-dir tmp -Oz -o drop_redundant_cnvs.sorted.vcf.gz drop_redundant_cnvs.vcf.gz
    tabix drop_redundant_cnvs.sorted.vcf.gz

    # --- StitchFragmentedCnvs ---
    java -Xmx${avail_mem}M -jar \$STITCH_JAR 0.2 200000 0.2 drop_redundant_cnvs.sorted.vcf.gz | bgzip > stitch1.vcf.gz
    java -Xmx${avail_mem}M -jar \$STITCH_JAR 0.2 200000 0.2 stitch1.vcf.gz | bgzip > stitched.vcf.gz

    # --- RescueMobileElementDeletions ---
    python - stitched.vcf.gz bnd_del.vcf.gz <<'PYEOF'
    import sys
    import pysam
    fin = pysam.VariantFile(sys.argv[1])
    fo = pysam.VariantFile(sys.argv[2], 'w', header=fin.header)
    for record in fin:
        if record.info['SVTYPE'] in ['BND'] and record.info['STRANDS'] == "+-" and record.chrom == record.info['CHR2'] and record.info['END2'] - record.start < 10000:
            record.info['SVLEN'] = record.info['END2'] - record.start
            fo.write(record)
    fin.close()
    fo.close()
    PYEOF
    tabix -p vcf bnd_del.vcf.gz
    svtk vcf2bed bnd_del.vcf.gz -i ALL --include-filters bnd_del.bed
    bgzip bnd_del.bed
    bedtools coverage -wo -a bnd_del.bed.gz -b ${line1_reference} | awk '{if (\$NF>.5) print}' | cut -f4 \\
        | sed -e 's/\$/\\tDEL\\tPASS\\toverlap_LINE1/' > manual_revise.MEI_DEL_from_BND.SVID_SVTYPE_FILTER_INFO.tsv
    bedtools coverage -wo -a bnd_del.bed.gz -b ${hervk_reference} | awk '{if (\$NF>.5) print}' | cut -f4 \\
        | sed -e 's/\$/\\tDEL\\tPASS\\toverlap_HERVK/' >> manual_revise.MEI_DEL_from_BND.SVID_SVTYPE_FILTER_INFO.tsv
    python - manual_revise.MEI_DEL_from_BND.SVID_SVTYPE_FILTER_INFO.tsv stitched.vcf.gz rescue_me_dels.vcf.gz <<'PYEOF'
    import sys
    import pysam
    reset = {}
    with open(sys.argv[1]) as f:
        for line in f:
            pin = line.strip().split()
            if pin[0] not in reset:
                reset[pin[0]] = pin[3]
    fin = pysam.VariantFile(sys.argv[2])
    fo = pysam.VariantFile(sys.argv[3], 'w', header=fin.header)
    for record in fin:
        if record.id in reset:
            del record.filter['UNRESOLVED']
            record.info['SVTYPE'] = 'DEL'
            record.info['SVLEN'] = record.info['END2'] - record.start
            record.stop = record.info['END2']
            record.info.pop("CHR2")
            record.info.pop("END2")
            record.info.pop("UNRESOLVED_TYPE")
            if reset[record.id] == 'overlap_LINE1':
                record.alts = ('<DEL:ME:LINE1>',)
            if reset[record.id] == 'overlap_HERVK':
                record.alts = ('<DEL:ME:HERVK>',)
        fo.write(record)
    fin.close()
    fo.close()
    PYEOF

    # --- AddHighFDRFilters ---
    python - rescue_me_dels.vcf.gz high_fdr_filtered.vcf.gz <<'PYEOF'
    import sys
    import pysam
    with pysam.VariantFile(sys.argv[1], 'r') as fin:
        header = fin.header
        header.add_line('##FILTER=<ID=HIGH_ALGORITHM_FDR,Description="Categories of variants with low precision including Wham-only deletions and certain Scramble SVAs">')
        with pysam.VariantFile(sys.argv[2], 'w', header=header) as fo:
            for record in fin:
                if (record.info['ALGORITHMS'] == ('wham',) and record.info['SVTYPE'] == 'DEL') or \\
                        (record.info['ALGORITHMS'] == ('scramble',) and record.info['HIGH_SR_BACKGROUND'] and record.alts == ('<INS:ME:SVA>',)):
                    record.filter.add('HIGH_ALGORITHM_FDR')
                fo.write(record)
    PYEOF

    # --- AddRetroDelFilters ---
    python ${retro_del_script} \\
        high_fdr_filtered.vcf.gz \\
        ${intron_reference} \\
        retro_del_filtered.vcf.gz \\
        --contig ${meta.contig}

    # --- FinalCleanup ---
    python /opt/sv-pipeline/04_variant_resolution/scripts/rename_after_vcfcluster.py \\
        --prefix ${id_prefix} \\
        retro_del_filtered.vcf.gz stdout \\
        | bcftools annotate --no-version -e 'SVTYPE=="CNV" && SVLEN<5000' -x INFO/MEMBERS -Oz -o final_cleanup.vcf.gz
    tabix final_cleanup.vcf.gz

    # --- FormatVcfToOutput ---
    python /opt/sv-pipeline/scripts/format_svtk_vcf_for_gatk.py \\
        --vcf final_cleanup.vcf.gz \\
        --out ${prefix}.vcf.gz \\
        --ploidy-table ${ploidy_table} \\
        --bothside-pass-list ${bothside_pass_list} \\
        --background-fail-list ${background_fail_list}
    tabix ${prefix}.vcf.gz
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo "" | gzip > ${prefix}.vcf.gz
    touch ${prefix}.vcf.gz.tbi
    """
}
