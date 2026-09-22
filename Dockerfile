# =============================================================================
# gamble-core Container Image for Routines / Cloud Workflows
# =============================================================================
FROM rocker/geospatial:4.4.0

LABEL maintainer="gamble-core"
LABEL description="Execution environment for GAMBLE prior land-use and livestock models"

# Prevent interactive prompts during apt install
ENV DEBIAN_FRONTEND=noninteractive

# System libraries for the C++ cores and linear algebra.
# libarmadillo-dev is deliberately NOT here: RcppArmadillo vendors its own Armadillo headers,
# so installing the system one adds a package and changes nothing.
RUN apt-get update && apt-get install -y --no-install-recommends \
    libopenblas-dev \
    liblapack-dev \
    ca-certificates lsb-release wget \
    && rm -rf /var/lib/apt/lists/*

# Arrow's C++ libraries come from APACHE'S OWN apt repository, not the Ubuntu archive.
# rocker/geospatial:4.4.0 is Ubuntu 22.04 (jammy), where `apt-get install libarrow-dev` fails
# with "Unable to locate package" -- which is exactly how this build died before the repository
# was added. arrow is not optional here: every driver reads the 1 km covariate parquet.
RUN apt-get update \
 && wget -q -O /tmp/apache-arrow-apt-source.deb \
      "https://packages.apache.org/artifactory/arrow/$(lsb_release --id --short | tr 'A-Z' 'a-z')/apache-arrow-apt-source-latest-$(lsb_release --codename --short).deb" \
 && apt-get install -y --no-install-recommends /tmp/apache-arrow-apt-source.deb \
 && apt-get update \
 && apt-get install -y --no-install-recommends libarrow-dev libparquet-dev \
 && rm -f /tmp/apache-arrow-apt-source.deb \
 && rm -rf /var/lib/apt/lists/*

# ONE dependency list for BOTH deployment paths: init.R is the contract and this just runs
# it, so the Accelerator predefined stack and this image can never drift apart. It is copied
# ALONE and run BEFORE the source, so editing a model file does not invalidate the package
# layer and trigger a full reinstall.
WORKDIR /app
COPY init.R /app/init.R
RUN Rscript /app/init.R

# Copy repository code
COPY . /app

# Ensure entrypoint is executable
RUN chmod +x /app/entrypoint.sh

# Default entrypoint for DeepOrigin / WKUBE routines
ENTRYPOINT ["./entrypoint.sh"]
