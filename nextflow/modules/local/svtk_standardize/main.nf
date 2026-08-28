//
// Standardize one caller's raw VCF into GATK-SV's common SV schema:
// rewrites the sample column to ped_id, restricts to primary contigs,
// applies minimum-size filtering, and sets standard INFO fields
// (SVTYPE/CHR2/END/STRANDS/SVLEN/ALGORITHMS). Wraps GATK-SV's own `svtk
// standardize` (src/svtk/), the actual mechanism their StandardizeVCFs
// task uses (wdl/PESRPreprocessing.wdl) -- not a script we can vendor
// like bin/*.py: svtk is a real Python package with a Cython extension
// and pybedtools/bedtools as a dependency, so this runs in GATK-SV's own
// published sv_pipeline_docker image rather than one we build ourselves.
//
// Replaces three separate point-fixes (modules/local/{manta,wham}/
// fix_sample_id's bcftools-reheader rename, and
// modules/local/vcf_primary_contigs_only's post-hoc bcftools filter) that
// were independently reconstructing pieces of what svtk standardize
// already does correctly in one step -- discovered via a real
// `Expected ALGORITHMS field` GATKException from SVCluster on real data,
// since none of those three fixes set ALGORITHMS at all. See
// wdl/PESRPreprocessing.wdl's StandardizeVCFs task, which this
// reproduces (svtk standardize, then bcftools sort + index).
//
process SVTK_STANDARDIZE {
    tag "$meta.id"
    label 'process_low'

    container 'us.gcr.io/broad-dsde-methods/gatk-sv/sv-pipeline:2026-06-02-v1.1-24f0fd08'

    input:
    tuple val(meta), path(vcf)
    path(contigs_fai)
    val(source)
    val(min_size)

    output:
    tuple val(meta), path("*.${source}.std.vcf.gz")    , emit: vcf
    tuple val(meta), path("*.${source}.std.vcf.gz.tbi"), emit: tbi

    when:
    task.ext.when == null || task.ext.when

    script:
    // source is included in the output filename (not just meta.id):
    // this process runs once per caller per sample, so the same sample
    // produces one *.std.vcf.gz per caller -- without source in the name,
    // Manta's and Wham's output for the same sample would collide (e.g.
    // under a shared publishDir, or an accidental filename clash locally).
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    svtk standardize \\
        --sample-names ${meta.ped_id} \\
        --prefix ${source}_${meta.ped_id} \\
        --contigs ${contigs_fai} \\
        --min-size ${min_size} \\
        ${vcf} tmp.vcf ${source}

    bcftools sort tmp.vcf -Oz -o ${prefix}.${source}.std.vcf.gz
    tabix -p vcf ${prefix}.${source}.std.vcf.gz
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo "" | gzip > ${prefix}.${source}.std.vcf.gz
    touch ${prefix}.${source}.std.vcf.gz.tbi
    """
}
