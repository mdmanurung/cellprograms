# CP-4 sharing/artifact sensitivity — frozen simulation harness.
# Drives the production program_sharing() on synthetic loading-level fits.
# See PREREGISTRATION.md for frozen config, scenarios, arms, decision rules.
#
# Usage:
#   Rscript run_sharing_sensitivity.R smoke  <out_dir>   # reps 1-3 per scenario
#   Rscript run_sharing_sensitivity.R confirm <out_dir>  # 20 reps per scenario
# CMP_MC_CORES sets parallelism (default 8).

suppressPackageStartupMessages(library(cellprograms))

## ---- universe ---------------------------------------------------------------

.universe <- function() {
  list(
    background = sprintf("G%03d", 1:300),
    ig = c(paste0("IGHV1-", 1:20), paste0("IGKV1-", 1:5), "IGKC",
           paste0("IGLC", 1:3), "JCHAIN", paste0("IGHV2-", 1:10)),
    tcr = c(paste0("TRBV", 1:10), paste0("TRAV", 1:10)),
    sex = c("XIST", "RPS4Y1", "DDX3Y", "UTY", "KDM5D", "EIF1AY",
            "NLGN4Y", "TMSB4Y", "USP9Y", "ZFY"),
    genotype = c("GSTM1", "GSTT1", "CCL3L1", "CCL4L1", "HLA-DRB3",
                 "HLA-DRB4", "HLA-DQA1", "HLA-DQB1", "ERAP2", "UTS2"),
    state = sprintf("ST%02d", 1:60)
  )
}
.universe_genes <- function(u) unlist(u, use.names = FALSE)

## ---- scenarios --------------------------------------------------------------

.scenarios <- function() {
  list(P1_shared_biology = list(shared = "state"),
       N0_unrelated      = list(shared = "none"),
       N1_ambient_ig     = list(shared = "ig"),
       N2_vdj_lineage    = list(shared = "tcr_lineage"),
       N3_sex_shared     = list(shared = "sex"),
       N4_genotype_shared = list(shared = "genotype"))
}

# Loading-level simulation: 4 types x 3 programs over the 400-gene universe.
# Support genes load ~ N(1.2, 0.1); off-support ~ N(0, 0.02). PSDs 0.08 on
# support, 0.15 elsewhere. Mono P3 is degenerate (~1e-18) by design.
sim_loadings <- function(shared, seed) {
  u <- .universe()
  genes <- .universe_genes(u)
  G <- length(genes)
  cts <- c("CD8T", "NK", "B", "Mono")
  set.seed(seed)
  shared_coef <- function(support) stats::setNames(rnorm(length(support), 1.2, 0.1), support)
  random_support <- function(n = 10) sample(u$background, n)
  W <- S <- vector("list", length(cts)); names(W) <- names(S) <- cts
  base_shared <- switch(shared,
    state = shared_coef(u$state), ig = shared_coef(u$ig),
    sex = shared_coef(u$sex), genotype = shared_coef(u$genotype),
    tcr_lineage = shared_coef(u$tcr), none = NULL)
  for (ct in cts) {
    M <- matrix(rnorm(G * 3, 0, 0.02), G, 3,
                dimnames = list(genes, paste0(ct, "_P", 1:3)))
    P <- matrix(0.15, G, 3, dimnames = dimnames(M))
    if (!is.null(base_shared) && (shared != "tcr_lineage" || ct == "CD8T")) {
      co <- base_shared + rnorm(length(base_shared), 0, 0.02)  # tiny type noise
      M[names(co), 1] <- co
      P[names(co), 1] <- 0.08
    } else {
      sup <- random_support()
      M[sup, 1] <- rnorm(length(sup), 1.2, 0.1)
      P[sup, 1] <- 0.08
    }
    for (k in 2:3) {  # independent background programs
      sup <- random_support()
      M[sup, k] <- rnorm(length(sup), 1.2, 0.1)
      P[sup, k] <- 0.08
    }
    if (ct == "Mono") {  # degenerate program fixture
      M[, 3] <- 1e-18
      P[, 3] <- 1e-16
    }
    W[[ct]] <- M; S[[ct]] <- P
  }
  list(W = W, S = S)
}

as_sharing_fit <- function(W, S) {
  cts <- names(W)
  structure(list(
    loadings = W,
    fits = stats::setNames(lapply(cts, function(ct) list(flash = list(F_psd = S[[ct]]))), cts),
    canonicalization = stats::setNames(lapply(cts, function(ct) list(scale_multiplier = rep(1, ncol(W[[ct]])))), cts),
    cell_types = cts),
    class = c("cell_program_fit", "cell_program_data", "list"))
}

## ---- arms -------------------------------------------------------------------

.empty_row <- function(status, message = "") {
  list(status = status, message = message, panel_n = NA_integer_, panel_md5 = "",
       n_degenerate = NA_integer_, n_matches = NA_integer_, n_sig_matches = NA_integer_,
       min_q = NA_real_, n_groups = 0L, group_size = NA_integer_, group_types = "",
       classification = "", p_all_shared = NA_real_, frac_top_shared = NA_real_,
       frac_top_private = NA_real_, artifact_frac_top = NA_real_, v_identity = NA,
       warnings = "")
}

# Run one arm of one replicate: panel subset -> production program_sharing().
# Returns list(rows = list of complete row lists, res = raw program_sharing result).
run_arm <- function(fit0, arm, seed, posthoc_drop = character(0)) {
  u <- .universe()
  universe <- .universe_genes(u)
  panel <- switch(arm,
    full = universe,
    prespecified = setdiff(universe, artifact_genes(universe, c("vdj", "genotype_marker", "sex"))),
    posthoc = setdiff(universe, posthoc_drop))
  fit <- fit0
  fit$loadings <- lapply(fit$loadings, function(M) M[panel, , drop = FALSE])
  fit$fits <- lapply(fit$fits, function(f) { f$flash$F_psd <- f$flash$F_psd[panel, , drop = FALSE]; f })
  warns <- character(0)
  res <- tryCatch(
    withCallingHandlers(
      program_sharing(fit, matches = NULL, alpha = 0.05, n_perm = 199L,
                      seed = seed, verbose = FALSE),
      warning = function(w) {
        warns <<- c(warns, conditionMessage(w))
        rst <- findRestart("mufflewarning")  # base warning()
        if (is.null(rst)) rst <- findRestart("muffle_warning")  # rlang/cli
        if (!is.null(rst)) invokeRestart(rst)
      }),
    error = function(e) e)
  degen <- cellprograms:::.degenerate_programs(fit)
  base <- list(panel_n = length(panel),
               panel_md5 = substr(as.character(.md5vec(panel)), 1, 10),
               n_degenerate = sum(lengths(degen)),
               warnings = paste(unique(warns), collapse = " | "))
  if (inherits(res, "error")) {
    row <- utils::modifyList(.empty_row("error", conditionMessage(res)), base)
    return(list(rows = list(row), res = NULL))
  }
  mt <- res$matches
  grp_rows <- list()
  if (!length(res$groups)) {
    grp_rows[[1]] <- utils::modifyList(.empty_row("ok"), base)
  } else {
    for (gi in seq_along(res$groups)) {
      g <- res$groups[[gi]]
      r <- res$results[[gi]]
      hits <- lapply(seq_len(nrow(g)), function(j) {
        Wv <- fit$loadings[[g$cell_type[j]]][, g$program_id[j]]
        top <- names(sort(abs(Wv), decreasing = TRUE))[1:12]
        mean(top %in% artifact_genes(panel, c("vdj", "genotype_marker", "sex")))
      })
      grp_rows[[gi]] <- utils::modifyList(.empty_row("ok"), c(base, list(
        n_groups = length(res$groups), group_size = nrow(g),
        group_types = paste(g$cell_type, collapse = ";"),
        classification = if (is.null(r)) "error" else r$classification,
        p_all_shared = if (is.null(r)) NA_real_ else r$p_all_shared,
        frac_top_shared = if (is.null(r)) NA_real_ else r$frac_top_shared,
        frac_top_private = if (is.null(r)) NA_real_ else r$frac_top_private,
        artifact_frac_top = mean(unlist(hits)),
        v_identity = if (is.null(r)) NA else isTRUE(all.equal(unname(r$V), diag(nrow(g)))))))
    }
  }
  match_fields <- list(n_matches = if (is.null(mt)) 0L else nrow(mt),
                       n_sig_matches = if (is.null(mt)) 0L else sum(mt$q < 0.05),
                       min_q = if (is.null(mt) || !nrow(mt)) NA_real_ else min(mt$q))
  grp_rows <- lapply(grp_rows, function(gr) utils::modifyList(gr, match_fields))
  list(rows = grp_rows, res = res)
}

.md5vec <- function(x) {
  f <- tempfile(); on.exit(unlink(f))
  writeLines(x, f); unname(tools::md5sum(f))
}

## ---- one replicate ----------------------------------------------------------

run_replicate <- function(scn_name, scn, rep_seed) {
  sim <- sim_loadings(scn$shared, rep_seed)
  fit0 <- as_sharing_fit(sim$W, sim$S)
  rows <- list()
  arm_res <- list()
  for (arm in c("full", "prespecified")) {
    ra <- run_arm(fit0, arm, rep_seed)
    arm_res[[arm]] <- ra$res
    for (x in ra$rows) {
      rows[[length(rows) + 1L]] <- c(list(scenario = scn_name, seed = rep_seed,
                                          arm = arm, posthoc_n_drop = NA_integer_), x)
    }
  }
  # posthoc arm: data-dependent exclusion driven by the FULL arm's flagged genes
  drop <- character(0)
  res_full <- arm_res$full
  if (!is.null(res_full) && length(res_full$groups)) {
    universe <- .universe_genes(.universe())
    for (g in res_full$groups) {
      for (j in seq_len(nrow(g))) {
        Wv <- fit0$loadings[[g$cell_type[j]]][, g$program_id[j]]
        top <- names(sort(abs(Wv), decreasing = TRUE))[1:12]
        drop <- union(drop, intersect(top, artifact_genes(universe, c("vdj", "genotype_marker", "sex"))))
      }
    }
  }
  ra <- run_arm(fit0, "posthoc", rep_seed, posthoc_drop = drop)
  for (x in ra$rows) {
    rows[[length(rows) + 1L]] <- c(list(scenario = scn_name, seed = rep_seed,
                                        arm = "posthoc", posthoc_n_drop = length(drop)), x)
  }
  do.call(rbind, lapply(rows, function(x) as.data.frame(x, stringsAsFactors = FALSE)))
}

## ---- Q1: two-axis distinction ----------------------------------------------

run_q1 <- function(out_dir) {
  cts <- c("CD8T", "NK", "B", "Mono")
  N <- 60L
  donors <- paste0("d", seq_len(N))
  mk <- function(correlated, seed) {
    set.seed(seed)
    sim <- sim_loadings("state", seed)
    z <- rnorm(N)
    scores <- lapply(cts, function(ct) {
      S <- matrix(rnorm(N * 3), N, 3, dimnames = list(donors, paste0(ct, "_P", 1:3)))
      S[, 1] <- if (correlated) z + rnorm(N, sd = 0.3) else rnorm(N)
      S
    })
    names(scores) <- cts
    fit <- as_sharing_fit(sim$W, sim$S)
    fit$scores <- scores
    fit$observation_ids <- donors
    fit$sample_metadata <- data.frame(observation_id = donors, stringsAsFactors = FALSE)
    fit
  }
  out <- list()
  for (regime in c("independent_scores", "correlated_scores")) {
    fit <- mk(regime == "correlated_scores", 77000)
    shr <- tryCatch(program_sharing(fit, matches = NULL, alpha = 0.05, n_perm = 199L,
                                    seed = 77000, verbose = FALSE), error = function(e) e)
    cmp <- tryCatch(find_coordinated_programs(fit, n_perm = 199L, min_shared = 10L,
                                              min_donors = 10L, ambient = "none",
                                              seed = 77000), error = function(e) e)
    out[[regime]] <- data.frame(
      regime = regime,
      sharing_status = if (inherits(shr, "error")) "error" else "ok",
      sharing_n_groups = if (inherits(shr, "error")) NA else length(shr$groups),
      sharing_classification = if (inherits(shr, "error") || !length(shr$groups)) "" else
        paste(vapply(shr$results, function(r) if (is.null(r)) "error" else r$classification, ""), collapse = ";"),
      cmp_status = if (inherits(cmp, "error")) "error" else "ok",
      cmp_n = if (inherits(cmp, "error") || is.null(cmp$cmps)) NA else length(unique(cmp$cmps$cmp)),
      stringsAsFactors = FALSE)
  }
  df <- do.call(rbind, out)
  write.csv(df, file.path(out_dir, "results_q1.csv"), row.names = FALSE)
  df
}

## ---- driver ------------------------------------------------------------------

main <- function(mode = c("smoke", "confirm"), out_dir) {
  mode <- match.arg(mode)
  n_reps <- if (mode == "smoke") 3L else 20L
  mc_cores <- as.integer(Sys.getenv("CMP_MC_CORES", "8"))
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  scns <- .scenarios()
  grid <- expand.grid(si = seq_along(scns), r = seq_len(n_reps))
  t_start <- Sys.time()
  message(sprintf("[%s] start %s | %d scenario-reps | cores=%d", mode,
                  format(t_start, tz = "UTC", usetz = TRUE), nrow(grid), mc_cores))
  res <- parallel::mclapply(seq_len(nrow(grid)), function(i) {
    si <- grid$si[i]; r <- grid$r[i]
    seed <- 31000L + si * 100L + r
    tryCatch(run_replicate(names(scns)[si], scns[[si]], seed),
             error = function(e) {
               row <- c(list(scenario = names(scns)[si], seed = seed, arm = NA_character_,
                             posthoc_n_drop = NA_integer_),
                        .empty_row("error", conditionMessage(e)))
               as.data.frame(row, stringsAsFactors = FALSE)
             })
  }, mc.cores = mc_cores, mc.preschedule = TRUE)
  df <- do.call(rbind, res)
  write.csv(df, file.path(out_dir, paste0("results_", mode, ".csv")), row.names = FALSE)
  q1 <- run_q1(out_dir)
  manifest <- list(mode = mode, n_reps = n_reps, alpha = 0.05, n_perm = 199,
                   scenarios = names(.scenarios()),
                   seed_rule = "31000 + scenario*100 + replicate",
                   mc_cores = mc_cores,
                   cellprograms_sha = system("git -C /workspace/repo rev-parse HEAD", intern = TRUE),
                   r_version = R.version.string,
                   started = format(t_start, tz = "UTC", usetz = TRUE),
                   finished = format(Sys.time(), tz = "UTC", usetz = TRUE),
                   wall_minutes = as.numeric(difftime(Sys.time(), t_start, units = "mins")))
  writeLines(jsonlite::toJSON(manifest, auto_unbox = TRUE, pretty = TRUE),
           file.path(out_dir, "run_manifest.json"))
  message(sprintf("[%s] wrote %d rows to %s in %.1f min", mode, nrow(df), out_dir,
                  manifest$wall_minutes))
  invisible(list(main = df, q1 = q1))
}

args <- commandArgs(trailingOnly = TRUE)
if (length(args) >= 2) main(args[1], args[2])
