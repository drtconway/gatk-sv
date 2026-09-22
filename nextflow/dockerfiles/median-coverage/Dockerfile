# Container for ../../bin/medianCoverage.R (adapted from GATK-SV's
# src/WGD/bin/medianCoverage.R -- see wdl/MedianCov.wdl's CalcMedCov task,
# which this reproduces, and bin/README.md for what was adapted and why).
#
# Third-party dependencies: CRAN's optparse (CLI parsing, from upstream),
# plus data.table and matrixStats (added for this pipeline's own
# fread()/colMedians()-based rewrite of the original read.table()/apply()
# computation -- see bin/README.md). Still base R + CRAN packages only,
# no bioconductor or compiled genomics library, so this keeps its own
# small image rather than pulling GATK-SV's sv_pipeline_qc_docker (a
# large, general-purpose image built for far more than this one script).
#
# Build and push:
#   docker build -t drtomc/gatk-sv-nf-median-coverage:0.2.0 .
#   docker push drtomc/gatk-sv-nf-median-coverage:0.2.0

FROM r-base:4.4.1

RUN Rscript -e 'install.packages(c("optparse", "data.table", "matrixStats"), repos="https://cloud.r-project.org")'
