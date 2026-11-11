FROM rocker/r-ver:4.4.1

ARG GIT_COMMIT
ENV DPYD_GIT_COMMIT=${GIT_COMMIT}

RUN apt-get update && apt-get install -y --no-install-recommends \
    jags git ca-certificates libcairo2-dev libfontconfig1-dev \
    libfreetype6-dev libpng-dev libjpeg-dev libtiff-dev \
    && rm -rf /var/lib/apt/lists/*

# The V4 execution record pins the statistical runtime. Install coda first so
# installing rjags cannot select a different dependency version.
RUN R -e "install.packages('remotes', repos='https://cloud.r-project.org/'); remotes::install_version('coda', version='0.19-4.1', repos='https://cloud.r-project.org/', upgrade='never'); remotes::install_version('rjags', version='4-17', repos='https://cloud.r-project.org/', upgrade='never'); install.packages(c('dplyr', 'ggplot2'), repos='https://cloud.r-project.org/'); remove.packages('remotes')"
RUN R -e "stopifnot(getRversion() == '4.4.1', as.character(rjags::jags.version()) == '4.3.0', packageVersion('rjags') == '4.17', packageVersion('coda') == '0.19.4.1', capabilities('cairo')); sessionInfo()"

WORKDIR /analysis
COPY R/ R/
COPY tests/ tests/
COPY data/ data/
COPY results/ results/
COPY README.md Dockerfile .dockerignore .gitignore ./
CMD ["Rscript", "R/run_scenario.R", "primary"]
