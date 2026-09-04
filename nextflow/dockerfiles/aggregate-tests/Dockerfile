# Container for the vendored GATK-SV script in ../../bin/aggregate.py
# (wraps GATK-SV's src/sv-pipeline/02_evidence_assessment/
# 02e_metric_aggregation/scripts/aggregate.py -- see
# wdl/GenerateBatchMetrics.wdl's AggregateTests task, which this
# reproduces).
#
# Separate from modules/local/format_svtk_vcf_for_gatk's
# sv-scripts image: this script needs numpy/pandas, which those three
# scripts don't -- same "small, purpose-built container per script"
# reasoning as dockerfiles/median-coverage, rather than growing one image
# to cover every vendored script's union of dependencies.
#
# Does NOT install svtk itself -- bin/aggregate.py inlines the one
# function (get_called_samples) it needed from svtk.utils, specifically to
# avoid pulling in pybedtools/bedtools for code this script never calls.
# See bin/README.md's own note on this before "fixing" that inlining.
#
# Build and push:
#   docker build -t drtomc/gatk-sv-nf-aggregate-tests:0.1.0 .
#   docker push drtomc/gatk-sv-nf-aggregate-tests:0.1.0

FROM python:3.12-slim

RUN pip install --no-cache-dir pysam==0.23.0 numpy pandas
