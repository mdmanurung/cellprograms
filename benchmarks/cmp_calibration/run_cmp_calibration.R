# CP-3A CMP calibration audit — frozen simulation harness.
# Drives the production find_coordinated_programs() on synthetic score-level
# fits. See PREREGISTRATION.md for frozen config, scenarios, decision rules.
#
# Usage:
#   Rscript run_cmp_calibration.R smoke  <out_dir>   # 20 reps, n_perm=199
#   Rscript run_cmp_calibration.R confirm <out_dir>  # frozen confirmation arm

suppressPackageStartupMessages(library(cellprograms))

## ---- generator -----------------------------------------------------------

# Simulate donor-level scores for C cell types x K programs over N donors.
# within_rho: within-cell-type equicorrelation (valid dependence).
# block_shared_rho: cross-type donor component shared within blocks (the
#   exchangeability violation: residualizing block MEANS does not remove it).
# planted: list of list(ct = <types>, rho = <strength>, signs = <+/-1>) —
#   a shared donor driver added to one program per listed cell type.
# blocks: named donor -> site vector. missingness: "none" or "shared_high".
sim_scores <- function(N, C, K, seed, within_rho = 0.3, blocks = NULL,
                       block_shared_rho = 0, planted = list(), planted_pairs = list(),
                       missingness = "none", zero_var = FALSE) {
  set.seed(seed)
  donors <- paste0("d", seq_len(N))
  if (is.null(blocks)) blocks <- stats::setNames(rep("all", N), donors)
  u <- matrix(rnorm(N * C), N, C, dimnames = list(donors, NULL))  # type-level donor factor
  if (block_shared_rho > 0) {
    b_eff <- rnorm(N)[match(blocks[donors], unique(blocks[donors]))]
    # shared across types within a block: same block deviate enters every type
    u <- sqrt(1 - block_shared_rho) * u + sqrt(block_shared_rho) *
      matrix(scale(b_eff)[, 1], N, C)
  }
  out <- vector("list", C)
  names(out) <- paste0("CT", seq_len(C))
  for (c in seq_len(C)) {
    E <- matrix(rnorm(N * K), N, K)
    Z <- sqrt(within_rho) * u[, c] + sqrt(1 - within_rho) * E
    colnames(Z) <- paste0("CT", c, "_P", seq_len(K))
    rownames(Z) <- donors
    out[[c]] <- Z
  }
  for (pl in planted) {
    w <- rnorm(N)
    for (i in seq_along(pl$ct)) {
      c <- pl$ct[i]
      s <- if (is.null(pl$signs)) 1 else pl$signs[i]
      out[[c]][, 1] <- sqrt(pl$rho) * (s * w) + sqrt(1 - pl$rho) * out[[c]][, 1]
    }
  }
  # pairwise drivers: independent latent per edge, enabling sign-inconsistent
  # triangles (a single shared driver is always sign-consistent by construction)
  for (pw in planted_pairs) {
    w <- rnorm(N)
    out[[pw$ct[1]]][, 1] <- sqrt(pw$rho) * (pw$signs[1] * w) + sqrt(1 - pw$rho) * out[[pw$ct[1]]][, 1]
    out[[pw$ct[2]]][, 1] <- sqrt(pw$rho) * (pw$signs[2] * w) + sqrt(1 - pw$rho) * out[[pw$ct[2]]][, 1]
  }
  if (zero_var && K >= 2) out[[1]][, 2] <- 1  # constant program
  if (identical(missingness, "shared_high")) {
    lat <- rowMeans(u)
    drop <- donors[lat >= quantile(lat, 0.75)]
    for (c in seq_len(C)) out[[c]][drop, ] <- NA_real_
  }
  out
}

# Minimal fit object carrying synthetic scores into the production caller.
as_cmp_fit <- function(scores, blocks) {
  ids <- unique(unlist(lapply(scores, rownames)))
  md <- data.frame(observation_id = ids,
                   site = unname(blocks[ids]), stringsAsFactors = FALSE)
  list(scores = scores,
       loadings = lapply(scores, function(S) matrix(0, 2, ncol(S))),  # non-NULL placeholder
       observation_ids = ids,
       sample_metadata = md,
       preprocessing = list(max_factors = 30L))
}

## ---- scenarios ------------------------------------------------------------

.scenarios <- function() {
  list(
    A0_sentinel  = list(design = "sentinel",  gen = list(within_rho = 0.3)),
    A0_realistic = list(design = "realistic", gen = list(within_rho = 0.3)),
    A1_sentinel  = list(design = "sentinel",  gen = list(within_rho = 0.3), blocks = "equal3"),
    A1_realistic = list(design = "realistic", gen = list(within_rho = 0.3), blocks = "equal3"),
    B1_realistic = list(design = "realistic", gen = list(within_rho = 0.3,
                        planted = list(list(ct = 1:3, rho = 0.6, signs = c(1, 1, 1))))),
    C1_realistic = list(design = "realistic", gen = list(within_rho = 0.3, block_shared_rho = 0.55),
                        blocks = "equal3"),
    C2_realistic = list(design = "realistic", gen = list(within_rho = 0.3), missingness = "shared_high"),
    C3_realistic = list(design = "realistic", gen = list(within_rho = 0.3), blocks = "unequal_singleton"),
    C4_realistic = list(design = "realistic", gen = list(within_rho = 0.3), zero_var = TRUE),
    D1_realistic = list(design = "realistic", gen = list(within_rho = 0.3,
                        planted = list(list(ct = 1:3, rho = 0.6, signs = c(1, 1, -1))))),
    D2_realistic = list(design = "realistic", gen = list(within_rho = 0.3,
                        planted_pairs = list(
                          list(ct = c(1, 2), rho = 0.6, signs = c(1, 1)),
                          list(ct = c(2, 3), rho = 0.6, signs = c(1, 1)),
                          list(ct = c(1, 3), rho = 0.6, signs = c(1, -1)))))
  )
}

.make_blocks <- function(kind, N, donors) {
  switch(kind,
    equal3 = stats::setNames(rep(c("s1", "s2", "s3"), length.out = N), donors),
    unequal_singleton = stats::setNames(
      c(rep("big", N - 9), rep("small", 8), "solo"), donors),
    stats::setNames(rep("all", N), donors)
  )
}

## ---- one replicate ----------------------------------------------------------

run_replicate <- function(scn_name, scn, rep_seed, n_perm, alpha = 0.05) {
  N <- if (scn$design == "sentinel") 30L else 100L
  C <- if (scn$design == "sentinel") 2L else 6L
  K <- if (scn$design == "sentinel") 4L else 8L
  donors <- paste0("d", seq_len(N))
  blocks <- .make_blocks(if (is.null(scn$blocks)) "none" else scn$blocks, N, donors)
  gen <- modifyList(list(N = N, C = C, K = K, seed = rep_seed, blocks = blocks),
                    scn$gen)
  scores <- do.call(sim_scores, gen)
  fit <- as_cmp_fit(scores, blocks)
  min_n <- if (scn$design == "sentinel") 10L else 20L
  t0 <- Sys.time()
  res <- tryCatch(
    find_coordinated_programs(fit, strata = "site", mode = "marginal",
                              ambient = "none", alpha = alpha, n_perm = n_perm,
                              tau = 2, link_threshold = 3,
                              min_shared = min_n, min_donors = min_n,
                              seed = rep_seed),
    error = function(e) e
  )
  elapsed <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  if (inherits(res, "error")) {
    return(list(scenario = scn_name, seed = rep_seed, status = "error",
                message = conditionMessage(res), elapsed = elapsed))
  }
  nodes <- res$nodes
  if (is.null(nodes) || !nrow(nodes)) {
    return(list(scenario = scn_name, seed = rep_seed, status = "ok", elapsed = elapsed,
                n_nodes = 0L, n_sig = 0L, any_sig = FALSE,
                null_n_sig = 0L, null_n = 0L, plant_n_sig = 0L, plant_n = 0L,
                n_cmps = 0L, cmp_membership = "", cmp_signs = ""))
  }
  node_key <- paste(nodes$cell_type, nodes$program_id, sep = "|")
  planted_keys <- character(0)
  for (pl in gen$planted %||% list()) {
    planted_keys <- c(planted_keys, paste0("CT", pl$ct, "|CT", pl$ct, "_P1"))
  }
  for (pw in gen$planted_pairs %||% list()) {
    planted_keys <- c(planted_keys,
                      paste0("CT", pw$ct[1], "|CT", pw$ct[1], "_P1"),
                      paste0("CT", pw$ct[2], "|CT", pw$ct[2], "_P1"))
  }
  planted_keys <- unique(planted_keys)
  is_planted <- node_key %in% planted_keys
  sig <- nodes$node_p < alpha
  cmps <- res$cmps
  list(scenario = scn_name, seed = rep_seed, status = "ok", elapsed = elapsed,
       n_nodes = nrow(nodes), n_sig = sum(sig), any_sig = any(sig),
       null_n_sig = sum(sig & !is_planted), null_n = sum(!is_planted),
       plant_n_sig = sum(sig & is_planted), plant_n = sum(is_planted),
       n_cmps = if (is.null(cmps)) 0L else length(unique(cmps$cmp)),
       cmp_membership = if (is.null(cmps) || !nrow(cmps)) "" else
         paste(paste0("CMP", cmps$cmp, ":", cmps$cell_type, "|", cmps$program_id), collapse = ";"),
       cmp_signs = if (is.null(cmps) || !nrow(cmps)) "" else
         paste(cmps$sign, collapse = ","))
}

`%||%` <- function(a, b) if (is.null(a)) b else a

## ---- driver ------------------------------------------------------------------

main <- function(mode = c("smoke", "confirm"), out_dir) {
  mode <- match.arg(mode)
  scns <- .scenarios()
  n_reps <- if (mode == "smoke") 20L else 500L
  n_perm <- if (mode == "smoke") 199L else 499L
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  results <- list()
  i <- 0L
  for (si in seq_along(scns)) {
    for (r in seq_len(n_reps)) {
      seed <- if (mode == "smoke") 9000L + r else 202610180L + si * 1000L + r
      i <- i + 1L
      results[[i]] <- run_replicate(names(scns)[si], scns[[si]], seed, n_perm)
      if (i %% 20 == 0) message(sprintf("[%s] %d/%d done", mode, i, length(scns) * n_reps))
    }
    # checkpoint per scenario
    saveRDS(results, file.path(out_dir, paste0("checkpoint_", mode, ".rds")))
  }
  df <- do.call(rbind, lapply(results, function(x) {
    if (x$status == "error") {
      return(data.frame(scenario = x$scenario, seed = x$seed, status = "error",
                        message = x$message, elapsed = x$elapsed,
                        n_nodes = NA, n_sig = NA, any_sig = NA,
                        null_n_sig = NA, null_n = NA, plant_n_sig = NA, plant_n = NA,
                        n_cmps = NA, cmp_membership = "", cmp_signs = "",
                        stringsAsFactors = FALSE))
    }
    data.frame(scenario = x$scenario, seed = x$seed, status = x$status,
               message = "", elapsed = x$elapsed, n_nodes = x$n_nodes,
               n_sig = x$n_sig, any_sig = x$any_sig,
               null_n_sig = x$null_n_sig, null_n = x$null_n,
               plant_n_sig = x$plant_n_sig, plant_n = x$plant_n,
               n_cmps = x$n_cmps, cmp_membership = x$cmp_membership,
               cmp_signs = x$cmp_signs, stringsAsFactors = FALSE)
  }))
  write.csv(df, file.path(out_dir, paste0("results_", mode, ".csv")), row.names = FALSE)
  manifest <- list(mode = mode, n_reps = n_reps, n_perm = n_perm, alpha = 0.05,
                   tau = 2, link_threshold = 3,
                   scenarios = names(scns),
                   seed_rule = if (mode == "smoke") "9000 + replicate" else "202610180 + scenario*1000 + replicate",
                   cellprograms_sha = system("git -C /workspace/repo rev-parse HEAD", intern = TRUE),
                   r_version = R.version.string,
                   started = format(Sys.time(), tz = "UTC", usetz = TRUE))
  writeLines(jsonlite::toJSON(manifest, auto_unbox = TRUE, pretty = TRUE),
             file.path(out_dir, "run_manifest.json"))
  message("wrote ", nrow(df), " replicate rows to ", out_dir)
  invisible(df)
}

args <- commandArgs(trailingOnly = TRUE)
if (length(args) >= 2) main(args[1], args[2])
