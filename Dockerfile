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

# Install required R packages:
# - Gibbs samplers & C++ bridge: Rcpp, RcppArmadillo, qs2, matrixStats, posterior, abind
# - Parallelism & progress: future, future.apply, progressr
# - Data & reporting: data.table, arrow, patchwork, base64enc, spdep
RUN R -e "install.packages(c( \
    'Rcpp', 'RcppArmadillo', 'data.table', 'qs2', 'matrixStats', \
    'posterior', 'abind', 'future', 'future.apply', 'progressr', \
    'patchwork', 'base64enc', 'spdep', 'arrow' \
), repos='https://cloud.r-project.org/', Ncpus=parallel::detectCores())"

WORKDIR /app

# Copy repository code
COPY . /app

# Ensure entrypoint is executable
RUN chmod +x /app/entrypoint.sh

# Default entrypoint for DeepOrigin / WKUBE routines
ENTRYPOINT ["./entrypoint.sh"]
