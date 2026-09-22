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

# MCRI's network intercepts outbound HTTPS via Zscaler -- without its
# root CA trusted, install.packages()'s HTTPS call to CRAN fails
# certificate verification, same underlying problem as ../../../Dockerfile's
# own npm ci (see its own note on this). r-base is Debian-based, and R's
# install.packages() uses the libcurl download method by default here
# (confirmed: capabilities("libcurl") is TRUE on this image), which
# already reads the system trust store update-ca-certificates
# populates -- so trusting the CA at the OS level should be enough on
# its own, unlike Node (which doesn't reliably pick up OS trust store
# updates and needs NODE_EXTRA_CA_CERTS pointed at it explicitly). Set
# CURL_CA_BUNDLE explicitly anyway, as a safety net in case some code
# path here doesn't fall back to the OS default the way plain libcurl
# does -- not confirmed to be necessary, just cheap insurance for a
# build step that's fully opaque until it either works or doesn't.
RUN apt-get update && apt-get install -y --no-install-recommends \
    ca-certificates \
    curl \
    && rm -rf /var/lib/apt/lists/* \
    && curl -fsSL -o /usr/local/share/ca-certificates/zscaler-root.crt \
    http://repo.mcri.edu.au/ZCC/Certs/ZscalerRootCerts/ZscalerRootCertificate-2048-SHA256.crt \
    && update-ca-certificates

ENV CURL_CA_BUNDLE=/etc/ssl/certs/ca-certificates.crt

RUN Rscript -e 'install.packages(c("optparse", "data.table", "matrixStats"), repos="https://cloud.r-project.org")'
