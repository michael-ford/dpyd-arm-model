#!/usr/bin/env Rscript
source("R/run_scenario.R")
publication_verify_inputs()
metadata <- list(analysis_commit = publication_git_commit(), docker_image = "runtime-recorded",
                 input_manifest_sha256 = sha256_file("data/input_checksums.csv"))
generate_core_checkpoint(read.csv("data/core/primary_raw.csv", stringsAsFactors = FALSE),
  "runs/checkpoints/core", metadata)
