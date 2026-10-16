## C2 registered check: optimal vs greedy program matching on stability refits.
## See PREREG.md for the registered protocol and decision rule.
##
## Usage:
##   Rscript run_check.R --fit <path/to/K30_fit.rds> --nboot 10 \
##                       --out <output_dir> [--seed-offset 0] [--frac 0.8]
##
## Resumable: each completed boot is checkpointed to <out>/boots/b<bb>.rds and
## skipped on rerun. Requires the `cellprograms` package (R CMD INSTALL the
## repo checkout first) and `clue` for the optimal assignment.
##
## Output per fit: <out>/discordance.rds with per-(boot, ct, program) rows and
## per-boot weight totals; summarize with summarize.R.

suppressPackageStartupMessages({
  library(cellprograms)
  library(clue)
})

args <- commandArgs(TRUE)
get_arg <- function(flag, default = NULL) {
  i <- match(paste0("--", flag), args)
  if (is.na(i)) return(default)
  args[i + 1L]
}
fit_path <- get_arg("fit"); stopifnot(!is.na(fit_path))
nboot <- as.integer(get_arg("nboot", "10"))
out_dir <- get_arg("out", "outputs")
frac <- as.numeric(get_arg("frac", "0.8"))
sim_threshold <- 0.5   # package default, per PREREG

dir.create(file.path(out_dir, "boots"), recursive = TRUE, showWarnings = FALSE)
fit <- readRDS(fit_path)
pre <- fit$preprocessing
dataset <- tools::file_path_sans_ext(basename(fit_path))

## Optimal assignment on |cosine|: partner vector for the ORIGINAL side
## (program_a index -> refit program index or NA). solve_LSAP assigns rows to
## columns with nrow <= ncol; transpose when the original side is larger.
.optimal_partners <- function(cos_sim) {
  S <- abs(cos_sim)
  partner <- rep(NA_integer_, nrow(S))
  if (nrow(S) == 0L || ncol(S) == 0L) return(partner)
  if (nrow(S) <= ncol(S)) {
    partner <- as.integer(clue::solve_LSAP(S, maximum = TRUE))   # row i -> column p(i)
  } else {
    a <- as.integer(clue::solve_LSAP(t(S), maximum = TRUE))      # refit j -> original a[j]
    partner[a] <- seq_along(a)
  }
  partner
}

.greedy_partners <- function(cos_sim) {
  gm <- cellprograms:::.greedy_match(cos_sim)
  partner <- rep(NA_integer_, nrow(cos_sim))
  if (nrow(gm)) partner[gm$a] <- gm$b
  partner
}

## post-threshold: partners below sim_threshold become NA
.threshold_partners <- function(partner, cos_sim, threshold) {
  has <- !is.na(partner)
  sim <- rep(0, length(partner))
  sim[has] <- abs(cos_sim[cbind(which(has), partner[has])])
  partner[has & sim < threshold] <- NA_integer_
  partner
}

run_boot <- function(b) {
  set.seed(b)  # donor subsamples; fit seed set inside as provenance$seed + b
  x_b <- fit[c("sample_metadata", "feature_metadata", "cell_types",
               "observation_ids", "shared_genes", "dims")]
  x_b$pseudobulk <- lapply(fit$pseudobulk, function(Y) {
    Y[sample.int(nrow(Y), max(5L, floor(frac * nrow(Y)))), , drop = FALSE]
  })
  class(x_b) <- "cell_program_data"
  fit_b <- tryCatch(
    fit_celltype_programs(
      x_b, loading_prior = pre$loading_prior, max_factors = pre$max_factors,
      center = pre$centered, scale = pre$scaled, features = pre$feature_sets,
      var_type = pre$var_type, backfit = pre$backfit, nullcheck = pre$nullcheck,
      covariates = pre$covariates, covariate_mode = pre$covariate_mode,
      flash_control = pre$flash_control, seed = fit$provenance$seed + b
    ) |> canonicalize_programs(),
    error = function(e) NULL
  )
  if (is.null(fit_b)) return(NULL)

  rows <- list(); w_greedy <- 0; w_opt <- 0
  for (ct in names(fit$fits)) {
    Wa <- fit$loadings[[ct]]; Wb <- fit_b$loadings[[ct]]
    K0 <- ncol(Wa)
    if (is.null(Wa) || is.null(Wb) || K0 == 0L || ncol(Wb) == 0L) next
    genes <- intersect(rownames(Wa), rownames(Wb))
    if (length(genes) < 10L) next
    Wa <- Wa[genes, , drop = FALSE]; Wb <- Wb[genes, , drop = FALSE]
    na <- sqrt(colSums(Wa^2)); nb <- sqrt(colSums(Wb^2))
    na[na == 0] <- 1; nb[nb == 0] <- 1
    cos_sim <- crossprod(Wa, Wb) / outer(na, nb)

    pg <- .threshold_partners(.greedy_partners(cos_sim), cos_sim, sim_threshold)
    po <- .threshold_partners(.optimal_partners(cos_sim), cos_sim, sim_threshold)
    w_greedy <- w_greedy + sum(abs(cos_sim[cbind(seq_along(pg), pg)]), na.rm = TRUE)
    w_opt <- w_opt + sum(abs(cos_sim[cbind(seq_along(po), po)]), na.rm = TRUE)

    changed <- ifelse(is.na(pg) & is.na(po), FALSE,
                      !is.na(pg) != !is.na(po) | (!is.na(pg) & !is.na(po) & pg != po))
    rows[[length(rows) + 1L]] <- data.frame(
      boot = b, cell_type = ct,
      program_a = colnames(Wa),
      partner_greedy = ifelse(is.na(pg), NA_character_, colnames(Wb)[pg]),
      partner_optimal = ifelse(is.na(po), NA_character_, colnames(Wb)[po]),
      discordant = changed, stringsAsFactors = FALSE)
  }
  do.call(rbind, rows)
}

for (b in seq_len(nboot)) {
  f_b <- file.path(out_dir, "boots", sprintf("b%03d.rds", b))
  if (file.exists(f_b)) { message("boot ", b, ": cached"); next }
  t0 <- Sys.time()
  res <- tryCatch(run_boot(b), error = function(e) {
    message("boot ", b, " failed: ", conditionMessage(e)); NULL
  })
  message(sprintf("boot %d/%d done in %.1f min", b, nboot,
                  as.numeric(difftime(Sys.time(), t0, units = "mins"))))
  if (!is.null(res)) saveRDS(res, f_b)
}

boots <- lapply(list.files(file.path(out_dir, "boots"), full.names = TRUE), readRDS)
all_rows <- do.call(rbind, boots)
saveRDS(all_rows, file.path(out_dir, "discordance.rds"))
rate <- mean(all_rows$discordant)
message(sprintf("[%s] pooled discordance: %.4f (%d/%d programs)",
                dataset, rate, sum(all_rows$discordant), nrow(all_rows)))
