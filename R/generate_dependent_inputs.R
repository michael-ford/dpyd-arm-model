#!/usr/bin/env Rscript
source("R/run_scenario.R")
publication_verify_inputs()
metadata <- list(analysis_commit = publication_git_commit(), docker_image = "runtime-recorded",
                 input_manifest_sha256 = sha256_file("data/input_checksums.csv"))
generate_dependent_checkpoint(getwd(), "runs/checkpoints/dependent", metadata,
  ld_second_target_policy = "both")
