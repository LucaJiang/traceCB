# R figures use the same exported configuration as the Python figures.
figure_path <- function(variable, kind = "file") {
  path <- Sys.getenv(variable, unset = "")
  guidance <- paste0(
    "Run 'source scripts/config.sh' in the terminal before starting R; check ",
    variable, " in scripts/config.sh."
  )
  if (!nzchar(path)) stop(guidance, call. = FALSE)
  path <- path.expand(path)
  if (kind == "output") {
    dir.create(path, recursive = TRUE, showWarnings = FALSE)
    if (!dir.exists(path)) stop(paste("Cannot create output directory:", path), call. = FALSE)
  } else if (kind == "directory") {
    if (!dir.exists(path)) stop(paste("Missing directory:", path, guidance), call. = FALSE)
  } else if (!file.exists(path) || dir.exists(path)) {
    stop(paste("Missing input file:", path, guidance), call. = FALSE)
  }
  normalizePath(path, mustWork = TRUE)
}

coloc_files <- function(directory) {
  files <- list.files(directory, pattern = "_coloc\\.csv$", full.names = TRUE)
  if (!length(files)) {
    stop(paste(
      "No *_coloc.csv files in", directory,
      "; check TRACECB_COLOC_DIR in scripts/config.sh or run bash scripts/run_colocalization.sh."
    ), call. = FALSE)
  }
  files
}
