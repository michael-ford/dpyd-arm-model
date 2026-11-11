# Public input integrity and execution provenance; no private checkout required.
publication_input_path <- function(output_root, path) {
  if (startsWith(path, "../../../data/")) {
    return(file.path(getwd(), sub("^\\.\\./\\.\\./\\.\\./", "", path)))
  }
  file.path(output_root, path)
}

publication_verify_inputs <- function(project_root = getwd()) {
  manifest <- read.csv(file.path(project_root, "data", "input_checksums.csv"), stringsAsFactors = FALSE)
  if (!nrow(manifest) || anyDuplicated(manifest$path)) stop("Invalid public input checksum manifest.", call. = FALSE)
  observed <- vapply(file.path(project_root, manifest$path), sha256_file, character(1))
  if (!identical(unname(observed), as.character(manifest$sha256))) {
    stop("Public input SHA-256 mismatch: ", paste(manifest$path[observed != manifest$sha256], collapse = ", "), call. = FALSE)
  }
  invisible(TRUE)
}

publication_run_root <- function(output_root) {
  root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
  reference <- file.path(root, "results", "analysis", c("core", "dependent"))
  if (normalizePath(output_root, winslash = "/", mustWork = TRUE) %in% reference) {
    file.path(root, "runs")
  } else file.path(output_root, "runs")
}

publication_git_commit <- function() {
  value <- if (nzchar(Sys.which("git"))) suppressWarnings(system2("git", c("rev-parse", "HEAD"), stdout = TRUE, stderr = FALSE)) else character()
  if (length(value) != 1L || !grepl("^[0-9a-f]{40}$", value)) value <- Sys.getenv("DPYD_GIT_COMMIT", "")
  if (length(value) != 1L || !grepl("^[0-9a-f]{40}$", value)) {
    stop("Git commit unavailable. Build with --build-arg GIT_COMMIT=\"$(git rev-parse HEAD)\".", call. = FALSE)
  }
  value
}

publication_write_metadata <- function(metadata, path) {
  metadata$r_version <- R.version.string
  metadata$r_platform <- R.version$platform
  metadata$jags_version <- as.character(rjags::jags.version())
  metadata$rjags_version <- as.character(utils::packageVersion("rjags"))
  metadata$coda_version <- as.character(utils::packageVersion("coda"))
  metadata$master_seed <- 12345L
  writeLines(capture.output(dput(R.version)), file.path(dirname(path), "R.version.txt"))
  write.csv(metadata, path, row.names = FALSE)
}

publication_relative_attempt <- function(path, output_root) {
  run_root <- publication_run_root(output_root)
  if (identical(run_root, file.path(getwd(), "runs"))) {
    substring(path, nchar(getwd()) + 2L)
  } else substring(path, nchar(output_root) + 2L)
}

# Tests copy the final CSVs into temporary checkpoints before mutating rows.
publication_copy_checkpoint_inputs <- function(fixture) {
  path <- file.path(fixture, "scenario_input_manifest.csv")
  manifest <- read.csv(path, stringsAsFactors = FALSE)
  dir.create(file.path(fixture, "inputs"), recursive = TRUE, showWarnings = FALSE)
  for (column in c("raw_input_path", "nma_input_path")) {
    for (index in seq_len(nrow(manifest))) {
      original <- publication_input_path(fixture, manifest[[column]][[index]])
      relative <- file.path("inputs", basename(original))
      if (!file.copy(original, file.path(fixture, relative))) stop("Cannot copy checkpoint fixture input.")
      manifest[[column]][[index]] <- relative
    }
  }
  write.csv(manifest, path, row.names = FALSE)
}
