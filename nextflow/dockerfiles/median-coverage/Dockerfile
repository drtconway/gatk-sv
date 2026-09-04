# Container for the vendored GATK-SV R script in ../../bin/medianCoverage.R
# (wraps GATK-SV's src/WGD/bin/medianCoverage.R -- see wdl/MedianCov.wdl's
# CalcMedCov task, which this reproduces).
#
# Only third-party dependency is CRAN's optparse -- base R otherwise, so
# this gets its own small image rather than pulling GATK-SV's
# sv_pipeline_qc_docker (a large, general-purpose image built for far more
# than this one script).
#
# Build and push:
#   docker build -t drtomc/gatk-sv-nf-median-coverage:0.1.0 .
#   docker push drtomc/gatk-sv-nf-median-coverage:0.1.0

FROM r-base:4.4.1

RUN Rscript -e 'install.packages("optparse", repos="https://cloud.r-project.org")'
