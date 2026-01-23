# R environment for DPYD Arm-Based Network Meta-Analysis
# Uses rocker/r-ver for reproducible R version

FROM rocker/r-ver:4.4.1

LABEL maintainer="NMA Analysis Team"
LABEL description="Arm-based NMA using pcnetmeta for DPYD variant analysis"

# Install system dependencies for JAGS and R packages
RUN apt-get update && apt-get install -y --no-install-recommends \
    jags \
    libcurl4-openssl-dev \
    libssl-dev \
    libxml2-dev \
    libfontconfig1-dev \
    libfreetype6-dev \
    libpng-dev \
    libtiff5-dev \
    libjpeg-dev \
    && rm -rf /var/lib/apt/lists/*

# -----------------------------------------------------------------------------
# Install R packages with pinned versions for reproducibility
# -----------------------------------------------------------------------------
# Why pin versions?
# - Ensures identical results across different builds and environments
# - Prevents breaking changes from package updates affecting analysis
# - Critical for scientific reproducibility of NMA results
# - Key packages (pcnetmeta, rjags, coda) directly affect statistical output
# -----------------------------------------------------------------------------

# Install remotes for version-specific installation
RUN R -e "install.packages('remotes', repos='https://cloud.r-project.org/')"

# Install core NMA packages with pinned versions
# These versions are current stable releases as of 2025-01
RUN R -e "remotes::install_version('rjags', version = '4-16', repos = 'https://cloud.r-project.org/')"
RUN R -e "remotes::install_version('coda', version = '0.19-4.1', repos = 'https://cloud.r-project.org/')"
RUN R -e "remotes::install_version('pcnetmeta', version = '2.8', repos = 'https://cloud.r-project.org/')"

# Install supporting packages (less critical for reproducibility, use latest)
RUN R -e "install.packages(c( \
    'readxl', \
    'dplyr', \
    'tidyr', \
    'ggplot2', \
    'knitr', \
    'rmarkdown' \
), repos='https://cloud.r-project.org/', Ncpus=4)"

# Install visualization packages for publication figures
RUN R -e "install.packages(c('igraph', 'ggraph', 'patchwork', 'scales', 'stringr'), repos='https://cloud.r-project.org/', Ncpus=4)"

# Verify JAGS connection
RUN R -e "library(rjags); cat('JAGS version:', .Call('get_version', PACKAGE='rjags'), '\n')"

# Set working directory
WORKDIR /analysis

# Copy analysis files
COPY R/ ./R/
COPY data/ ./data/

# Default command: run the analysis
CMD ["Rscript", "R/run_analysis.R"]
