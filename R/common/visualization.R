# =============================================================================
# R/common/visualization.R
# Shared visualization utilities for NMA results
# =============================================================================

library(pcnetmeta)

# -----------------------------------------------------------------------------
# Generate all standard NMA visualizations
# -----------------------------------------------------------------------------

generate_visualizations <- function(result, nma_data, treatment_names,
                                     model_name, output_dir = "output") {
  cat("\n=== Generating NMA Visualizations ===\n\n")

  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }

  reference_name <- treatment_names[1]

  # 1. Network Plot
  cat("Creating network plot...\n")
  tryCatch({
    filename <- file.path(output_dir, paste0(model_name, "_network_plot.png"))
    png(filename, width = 800, height = 800, res = 120)
    nma.networkplot(
      s.id = s.id,
      t.id = t.id,
      data = nma_data,
      trtname = treatment_names,
      title = paste("Network Structure:", model_name)
    )
    dev.off()
    cat("  Saved:", filename, "\n")
  }, error = function(e) {
    cat("  Error creating network plot:", e$message, "\n")
    try(dev.off(), silent = TRUE)
  })

  # 2. Contrast Plot (ORs vs reference)
  cat("Creating contrast plot...\n")
  tryCatch({
    filename <- file.path(output_dir, paste0(model_name, "_contrast_plot.png"))
    png(filename, width = 1000, height = 600, res = 120)
    contrast.plot(
      result,
      reference = reference_name,
      digits = 2
    )
    dev.off()
    cat("  Saved:", filename, "\n")
  }, error = function(e) {
    cat("  Error creating contrast plot:", e$message, "\n")
    try(dev.off(), silent = TRUE)
  })

  # 3. Absolute Effects Plot
  cat("Creating absolute effects plot...\n")
  tryCatch({
    filename <- file.path(output_dir, paste0(model_name, "_absolute_effects.png"))
    png(filename, width = 1000, height = 600, res = 120)
    absolute.plot(result, digits = 2)
    dev.off()
    cat("  Saved:", filename, "\n")
  }, error = function(e) {
    cat("  Error creating absolute effects plot:", e$message, "\n")
    try(dev.off(), silent = TRUE)
  })

  # 4. Ranking Plot
  cat("Creating ranking plot...\n")
  tryCatch({
    filename <- file.path(output_dir, paste0(model_name, "_ranking_plot.png"))
    png(filename, width = 1000, height = 600, res = 120)
    rank.prob(result)
    dev.off()
    cat("  Saved:", filename, "\n")
  }, error = function(e) {
    cat("  Error creating ranking plot:", e$message, "\n")
    try(dev.off(), silent = TRUE)
  })

  cat("\n=== Visualizations Complete ===\n")
  cat("Generated plots:\n")
  cat("  -", paste0(model_name, "_network_plot.png"), "(network structure)\n")
  cat("  -", paste0(model_name, "_contrast_plot.png"), "(ORs vs", reference_name, ")\n")
  cat("  -", paste0(model_name, "_absolute_effects.png"), "(toxicity rates)\n")
  cat("  -", paste0(model_name, "_ranking_plot.png"), "(treatment rankings)\n")
}

# -----------------------------------------------------------------------------
# Move trace plots to output directory
# -----------------------------------------------------------------------------

organize_trace_plots <- function(model_name, output_dir = "output") {
  # pcnetmeta generates trace plots in working directory
  # Move them to the output directory with model prefix

  trace_files <- list.files(pattern = "^LOR.*\\.png$")
  if (length(trace_files) > 0) {
    cat("\nOrganizing trace plots...\n")
    for (f in trace_files) {
      new_name <- paste0(model_name, "_trace_", f)
      new_path <- file.path(output_dir, new_name)
      file.rename(f, new_path)
      cat("  Moved:", f, "->", new_name, "\n")
    }
  }

  # Also move convergence diagnostic file
  conv_files <- list.files(pattern = "ConvergenceDiagnostic.*\\.txt$")
  if (length(conv_files) > 0) {
    for (f in conv_files) {
      new_name <- paste0(model_name, "_", f)
      new_path <- file.path(output_dir, new_name)
      file.rename(f, new_path)
      cat("  Moved:", f, "->", new_name, "\n")
    }
  }
}
