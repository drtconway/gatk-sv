# Container for the vendored svtk adjudicate logic in ../../bin/
# (adjudicate.py, adjudicate_sv.py, random_forest.py, labelers.py --
# wraps GATK-SV's svtk.adjudicate package, see wdl/FilterBatchSites.wdl's
# AdjudicateSV task, which this reproduces).
#
# Only third-party deps are pandas/numpy/scikit-learn -- no pysam/
# pybedtools/bedtools/Cython needed, unlike a real `svtk` install: the
# actual svtk.adjudicate code has no such dependency, but invoking it via
# svtk's own CLI dispatcher (scripts/svtk) would need the whole package
# installed anyway, since that dispatcher eagerly imports every
# subcommand module at startup (several of which do need pysam/
# pybedtools). Vendored and called directly instead, to avoid that --
# see bin/README.md's own note.
#
# Build and push:
#   docker build -t drtomc/gatk-sv-nf-adjudicate:0.1.0 .
#   docker push drtomc/gatk-sv-nf-adjudicate:0.1.0

FROM python:3.12-slim

RUN pip install --no-cache-dir pandas numpy scikit-learn
