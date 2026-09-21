# =============================================================================
# gamble-core Container Image for Routines / Cloud Workflows
# =============================================================================
FROM rocker/geospatial:4.4.0

LABEL maintainer="gamble-core"
LABEL description="Execution environment for GAMBLE prior land-use and livestock models"

# Prevent interactive prompts during apt install
ENV DEBIAN_FRONTEND=noninteractive

# Install system libraries needed for C++ cores, parquet/arrow, and linear algebra
RUN apt-get update && apt-get install -y --no-install-recommends \
    libarmadillo-dev \
    libarrow-dev \
    libparquet-dev \
    libopenblas-dev \
    liblapack-dev \
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
