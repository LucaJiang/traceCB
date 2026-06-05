#!/usr/bin/env Rscript

print_usage <- function() {
  cat(
    "Usage: Rscript src/simulation/others/run_mashr_benchmark.R [options]\n\n",
    "Options:\n",
    "  --base_path PATH          Base directory containing simulation outputs\n",
    "  --runname NAME            Simulation run name under base_path\n",
    "  --cov_method METHOD       mashr covariance method\n",
    "  --lfsr_threshold VALUE    lfsr significance threshold\n",
    "  --pvalue_threshold VALUE  Posterior-mean z-test p-value threshold\n",
    "  --pca_factors N           Number of PCA factors for canonical_pca\n",
    "  --condition_set SET       One of sc or sc_bulk\n",
    "  --output_prefix PREFIX    Output CSV prefix written beside each simulation file\n",
    "  --max_files N             Limit number of simulation files processed\n",
    "  --force                   Recompute outputs even if they already exist\n",
    "  -h, --help                Show this help message\n",
    sep = ""
  )
}

parse_args <- function(args) {
  opts <- list(
    base_path = "bench/result_estOmega",
    runname = "alpha_h2sq_pcausal_propt",
    cov_method = "canonical_pca",
    lfsr_threshold = 0.05,
    pvalue_threshold = 0.05,
    pca_factors = 3,
    condition_set = "sc_bulk",
    output_prefix = "mashr",
    max_files = Inf,
    force = FALSE
  )
  i <- 1
  while (i <= length(args)) {
    key <- args[[i]]
    if (key %in% c("--help", "-h")) {
      print_usage()
      quit(save = "no", status = 0)
    }
    if (key == "--force") {
      opts$force <- TRUE
      i <- i + 1
      next
    }
    if (i == length(args)) {
      stop(paste("Missing value for", key))
    }
    value <- args[[i + 1]]
    if (key == "--base_path") opts$base_path <- value
    else if (key == "--runname") opts$runname <- value
    else if (key == "--cov_method") opts$cov_method <- value
    else if (key == "--lfsr_threshold") opts$lfsr_threshold <- as.numeric(value)
    else if (key == "--pvalue_threshold") opts$pvalue_threshold <- as.numeric(value)
    else if (key == "--pca_factors") opts$pca_factors <- as.integer(value)
    else if (key == "--condition_set") opts$condition_set <- value
    else if (key == "--output_prefix") opts$output_prefix <- value
    else if (key == "--max_files") opts$max_files <- as.integer(value)
    else stop(paste("Unknown option", key))
    i <- i + 2
  }
  opts
}

parse_setting <- function(setting_dir) {
  tokens <- strsplit(basename(setting_dir), "_", fixed = TRUE)[[1]]
  params <- list()
  i <- 1
  while (i <= length(tokens)) {
    key <- tokens[[i]]
    if (key %in% c("h1sq", "h2sq", "gc", "n1", "n2", "nt", "nsnp", "propt", "pcausal", "omega") &&
        i + 1 <= length(tokens)) {
      raw <- tokens[[i + 1]]
      if (tolower(raw) == "true") {
        value <- TRUE
      } else if (tolower(raw) == "false") {
        value <- FALSE
      } else {
        value <- as.numeric(raw)
      }
      params[[key]] <- value
      i <- i + 2
    } else {
      i <- i + 1
    }
  }
  for (required in c("n1", "n2", "nt")) {
    if (is.null(params[[required]])) {
      stop(paste("Could not parse", required, "from", basename(setting_dir)))
    }
  }
  params
}

symmetrize_Ulist <- function(Ulist) {
  lapply(Ulist, function(U) {
    U <- as.matrix(U)
    (U + t(U)) / 2
  })
}

run_mash_with_fallback <- function(mash_data, Ulist) {
  attempts <- list(
    default = list(),
    mixEM = list(optmethod = "mixEM"),
    R_version = list(algorithm.version = "R")
  )
  errors <- c()
  for (attempt_name in names(attempts)) {
    args <- c(
      list(data = mash_data, Ulist = Ulist, verbose = FALSE),
      attempts[[attempt_name]]
    )
    fit <- tryCatch(
      do.call(mash, args),
      error = function(e) e
    )
    if (!inherits(fit, "error")) {
      return(list(fit = fit, status = attempt_name, error = ""))
    }
    errors <- c(errors, paste0(attempt_name, ": ", conditionMessage(fit)))
  }
  list(fit = NULL, status = "failed", error = paste(errors, collapse = " | "))
}

get_condition_spec <- function(condition_set, setting) {
  if (condition_set == "sc") {
    return(list(
      conditions = c("pop1sc", "pop2sc"),
      z_cols = c("z1_sumstat", "z2_sumstat"),
      b_cols = c("b1_hat", "b2_hat"),
      se_cols = c("se1_hat", "se2_hat"),
      se = c(
        pop1sc = 1 / sqrt(as.numeric(setting$n1)),
        pop2sc = 1 / sqrt(as.numeric(setting$n2))
      )
    ))
  }
  if (condition_set == "sc_bulk") {
    return(list(
      conditions = c("pop1sc", "pop2sc", "pop2bulk_tissue"),
      z_cols = c("z1_sumstat", "z2_sumstat", "zt_sumstat"),
      b_cols = c("b1_hat", "b2_hat", "bt_hat"),
      se_cols = c("se1_hat", "se2_hat", "se_t_hat"),
      se = c(
        pop1sc = 1 / sqrt(as.numeric(setting$n1)),
        pop2sc = 1 / sqrt(as.numeric(setting$n2)),
        pop2bulk_tissue = 1 / sqrt(as.numeric(setting$nt))
      )
    ))
  }
  stop(paste("Unknown condition_set", condition_set))
}

make_mash_data_from_simulation <- function(simulation_file, condition_set) {
  setting <- parse_setting(dirname(simulation_file))
  dat <- read.csv(
    simulation_file,
    check.names = FALSE,
    colClasses = c(causal = "numeric")
  )
  spec <- get_condition_spec(condition_set, setting)
  if (all(c(spec$b_cols, spec$se_cols) %in% names(dat))) {
    Bhat <- as.matrix(dat[, spec$b_cols])
    Shat <- as.matrix(dat[, spec$se_cols])
  } else {
    z <- as.matrix(dat[, spec$z_cols])
    Bhat <- sweep(z, 2, spec$se, "*")
    Shat <- matrix(rep(spec$se, each = nrow(z)), nrow = nrow(z), byrow = FALSE)
  }
  colnames(Bhat) <- spec$conditions
  colnames(Shat) <- spec$conditions
  list(
    Bhat = Bhat,
    Shat = Shat,
    causal = dat$causal,
    conditions = spec$conditions
  )
}

fit_one_file <- function(simulation_file, opts) {
  output_file <- sub(
    "simulation_([0-9]+)\\.csv$",
    paste0(opts$output_prefix, "_\\1.csv"),
    simulation_file
  )
  if (!opts$force && file.exists(output_file)) {
    return(FALSE)
  }

  sim <- make_mash_data_from_simulation(simulation_file, opts$condition_set)
  mash_data <- mash_set_data(sim$Bhat, sim$Shat)
  Ulist <- cov_canonical(mash_data)

  strong_n <- 0
  if (opts$cov_method %in% c("canonical_pca", "canonical_pca_ed")) {
    m1by1 <- tryCatch(
      mash_1by1(mash_data),
      error = function(e) {
        message("mash_1by1 failed for ", simulation_file, ": ", conditionMessage(e))
        NULL
      }
    )
    if (!is.null(m1by1)) {
      strong <- get_significant_results(m1by1, opts$lfsr_threshold)
      strong_n <- length(strong)
      if (strong_n >= 2) {
        npcs <- min(opts$pca_factors, ncol(sim$Bhat), strong_n)
        U_pca <- tryCatch(
          cov_pca(mash_data, npcs, subset = strong),
          error = function(e) NULL
        )
        if (!is.null(U_pca)) {
          Ulist <- c(Ulist, U_pca)
          if (opts$cov_method == "canonical_pca_ed") {
            U_ed <- tryCatch(
              cov_ed(mash_data, U_pca, subset = strong),
              error = function(e) NULL
            )
            if (!is.null(U_ed)) Ulist <- c(Ulist, U_ed)
          }
        }
      }
    }
  } else if (opts$cov_method != "canonical") {
    stop(paste("Unknown cov_method", opts$cov_method))
  }

  Ulist <- symmetrize_Ulist(Ulist)
  fit_result <- run_mash_with_fallback(mash_data, Ulist)
  if (is.null(fit_result$fit)) {
    message("mash failed for ", simulation_file, ": ", fit_result$error)
    out <- data.frame(
      snp_id = seq_len(nrow(sim$Bhat)) - 1,
      causal = sim$causal,
      condition_set = opts$condition_set,
      mash_lfsr_pop1sc = 1.0,
      z1_mashr_pm = NA_real_,
      p1_mashr_pm = 1.0,
      sig1_mashr_lfsr = FALSE,
      sig1_mashr_pm = FALSE,
      mash_strong_n = strong_n,
      mash_fit_status = "failed",
      mash_error = fit_result$error
    )
    for (condition in setdiff(sim$conditions, "pop1sc")) {
      out[[paste0("mash_lfsr_", condition)]] <- 1.0
      out[[paste0("z_mashr_pm_", condition)]] <- NA_real_
      out[[paste0("p_mashr_pm_", condition)]] <- 1.0
    }
    if ("pop2sc" %in% sim$conditions) {
      out$z2_mashr_pm <- NA_real_
      out$p2_mashr_pm <- 1.0
    }
    if ("pop2bulk_tissue" %in% sim$conditions) {
      out$zt_mashr_pm <- NA_real_
      out$pt_mashr_pm <- 1.0
    }
    write.csv(out, output_file, row.names = FALSE)
    return(TRUE)
  }
  fit <- fit_result$fit
  lfsr <- get_lfsr(fit)
  pm <- get_pm(fit)
  psd <- get_psd(fit)
  z_pm <- pm / pmax(psd, .Machine$double.eps)
  p_pm <- 2 * pnorm(abs(z_pm), lower.tail = FALSE)

  out <- data.frame(
    snp_id = seq_len(nrow(sim$Bhat)) - 1,
    causal = sim$causal,
    condition_set = opts$condition_set,
    mash_lfsr_pop1sc = lfsr[, "pop1sc"],
    z1_mashr_pm = z_pm[, "pop1sc"],
    p1_mashr_pm = p_pm[, "pop1sc"],
    sig1_mashr_lfsr = lfsr[, "pop1sc"] < opts$lfsr_threshold,
    sig1_mashr_pm = p_pm[, "pop1sc"] < opts$pvalue_threshold,
    mash_strong_n = strong_n,
    mash_fit_status = fit_result$status,
    mash_error = fit_result$error
  )
  for (condition in setdiff(sim$conditions, "pop1sc")) {
    out[[paste0("mash_lfsr_", condition)]] <- lfsr[, condition]
    out[[paste0("z_mashr_pm_", condition)]] <- z_pm[, condition]
    out[[paste0("p_mashr_pm_", condition)]] <- p_pm[, condition]
  }
  if ("pop2sc" %in% sim$conditions) {
    out$z2_mashr_pm <- z_pm[, "pop2sc"]
    out$p2_mashr_pm <- p_pm[, "pop2sc"]
  }
  if ("pop2bulk_tissue" %in% sim$conditions) {
    out$zt_mashr_pm <- z_pm[, "pop2bulk_tissue"]
    out$pt_mashr_pm <- p_pm[, "pop2bulk_tissue"]
  }
  write.csv(out, output_file, row.names = FALSE)
  TRUE
}

opts <- parse_args(commandArgs(trailingOnly = TRUE))
suppressPackageStartupMessages(library(mashr))
run_dir <- file.path(opts$base_path, opts$runname)
if (!dir.exists(run_dir)) {
  stop(paste("Run directory does not exist:", run_dir))
}
files <- sort(list.files(
  run_dir,
  pattern = "^simulation_[0-9]+\\.csv$",
  recursive = TRUE,
  full.names = TRUE
))
if (!is.infinite(opts$max_files)) {
  files <- head(files, opts$max_files)
}
if (length(files) == 0) {
  stop(paste("No simulation files found under", run_dir))
}

cat("mashr benchmark start:", date(), "\n")
cat(
  "Simulation files:", length(files),
  "cov_method:", opts$cov_method,
  "condition_set:", opts$condition_set,
  "output_prefix:", opts$output_prefix,
  "\n"
)
written <- 0
for (i in seq_along(files)) {
  wrote <- fit_one_file(files[[i]], opts)
  written <- written + as.integer(wrote)
  if (i %% 20 == 0 || i == length(files)) {
    cat("mashr files", i, "/", length(files), "done; written", written, "\n")
  }
}
cat("mashr benchmark end:", date(), "\n")
