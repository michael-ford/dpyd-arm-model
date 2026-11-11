#!/usr/bin/env Rscript
# Rebuild publication assets from committed lightweight numerical results.
source("R/checksums.R")
publication_verify_inputs()
source("R/render_core_figures.R")
source("R/build_core_tables.R")
source("R/render_sensitivity_figures.R")
source("R/build_sensitivity_tables.R")
render_core_figures()
build_core_tables()
render_sensitivity_figures()
build_sensitivity_tables()
cat("Rebuilt tables and figures in runs/rebuilt/.\n")
