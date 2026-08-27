# Container for the vendored GATK-SV Python scripts in ../../bin/
# (format_svtk_vcf_for_gatk.py, format_gatk_vcf_for_svtk.py,
# ploidy_table_from_ped.py) that need pysam + tabix.
#
# Known gap: uses current pysam, not GATK-SV's pinned pysam==0.15.4 --
# see ../../bin/README.md for the BND/CTX record risk this may carry.
#
# Build and push:
#   docker build -t drtomc/gatk-sv-nf-sv-scripts:0.1.0 .
#   docker push drtomc/gatk-sv-nf-sv-scripts:0.1.0
#
# Referenced by modules/local/format_svtk_vcf_for_gatk/main.nf and
# modules/local/format_gatk_vcf_for_svtk/main.nf.

FROM python:3.12-slim

RUN apt-get update -qq \
    && apt-get install -qq -y --no-install-recommends tabix \
    && rm -rf /var/lib/apt/lists/*

RUN pip install --no-cache-dir pysam==0.23.0
