# Container for FILTER_ANNOTATE_VCF (modules/local/filter_annotate_vcf),
# wrapping GATK-SV's FilterBatchSites.wdl FilterAnnotateVcf task. Needs
# both bcftools/bgzip (the task's own inline shell filtering pipeline)
# and pysam/numpy/pandas (../../bin/rewrite_SR_coords.py,
# ../../bin/annotate_RF_evidence.py) in the same container, since a
# single Nextflow process has exactly one container -- neither this
# pipeline's existing bcftools_htslib image nor its aggregate-tests image
# (python:3.12-slim + pysam/numpy/pandas, no bcftools) covers both, so
# this is a new, slightly larger combined image rather than splitting one
# WDL task's single inline pipeline across two Nextflow processes.
#
# Build and push:
#   docker build -t drtomc/gatk-sv-nf-filter-annotate-vcf:0.1.0 .
#   docker push drtomc/gatk-sv-nf-filter-annotate-vcf:0.1.0

FROM python:3.12-slim

RUN apt-get update -qq \
    && apt-get install -qq -y --no-install-recommends bcftools tabix \
    && rm -rf /var/lib/apt/lists/*

RUN pip install --no-cache-dir pysam==0.23.0 numpy pandas
