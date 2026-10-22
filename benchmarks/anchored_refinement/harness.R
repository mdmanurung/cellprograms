## harness.R -- MR-4 anchored refinement + ablation (REG-MR4-anchored-v1).
## Anchor = Stage A v0.2.0 fits (donor pseudobulk, K30 recipe); query =
## cell-level (or metacell-level) counts; anchored arms refit with
## flashier.utils::fit_initialized_ebmf(). Sources the repo's R/*.R.
##
## Conventions: X_cell = cells x genes (samples in rows, sample_side = "rows");
## activities = cells x factors; loadings = genes x factors.

library(Matrix)

`%||%` <- function(a, b) if (is.null(a)) b else a

## repo root: CP_REPO env, else two levels up from the working directory
REPO <- Sys.getenv("CP_REPO", unset = normalizePath(file.path(getwd(), "..", ".."), mustWork = FALSE))
stopifnot(dir.exists(file.path(REPO, "R")))
for (f in sort(list.files(file.path(REPO, "R"), pattern = "\\.R$", full.names = TRUE))) source(f)

## ---- simulator: semi-synthetic cell-level counts ---------------------------
## 4 cell types (A-D), n_donors donors (2 sites), G genes. A 50-gene shared
## program is planted in the A-B pair (donor activities z ~ N(0,1)); 2 private
## programs per type. Cells per donor lognormal(median, sd_log); Poisson
## counts on exp(log-profile). Returns cell-level counts + truth.
sim_cell_level <- function(n_donors = 60, G = 5000, program_size = 50,
                           amp = 1, n_cells_med = 200, seed = 1, n_ct = 4L) {
  set.seed(seed)
  donors <- sprintf("d%03d", seq_len(n_donors))
  site <- setNames(rep(c("s1", "s2"), length.out = n_donors), donors)
  cts <- c("A", "B", "C", "D")[seq_len(n_ct)]
  base_mu <- 2 + 2 * rexp(G, rate = 1)                 # gene means (log scale)
  site_eff <- matrix(rnorm(G * 2, sd = 0.15), G, 2, dimnames = list(NULL, c("s1", "s2")))
  donor_eff <- matrix(rnorm(G * n_donors, sd = 0.2), G, n_donors, dimnames = list(NULL, donors))
  shared_genes <- sort(sample(G, program_size))
  shared_z <- setNames(rnorm(n_donors), donors)        # the planted CMP gradient
  priv_genes <- stats::setNames(lapply(cts, function(ct) setdiff(sample(G, program_size * 2), shared_genes)), cts)
  priv_z <- stats::setNames(lapply(cts, function(ct) setNames(rnorm(n_donors), donors)), cts)
  ind_shared <- numeric(G); ind_shared[shared_genes] <- 1
  ind_priv <- stats::setNames(lapply(priv_genes, function(g) { v <- numeric(G); v[g] <- 1; v }), cts)
  rows <- list(); meta <- list()
  n_i <- setNames(pmax(20, round(exp(rnorm(n_donors, log(n_cells_med), 0.8)))), donors)
  for (ct in cts) {
    for (d in donors) {
      logprof <- rep(base_mu + site_eff[, site[d]] + donor_eff[, d] +
                       amp * shared_z[d] * ind_shared * (ct %in% c("A", "B")) +
                       amp * priv_z[[ct]][d] * ind_priv[[ct]], n_i[d])
      ## one independent Poisson draw per gene per cell (cell-blocked lambda);
      ## clip at a biologically plausible max count (heavy-tailed base_mu can
      ## otherwise overflow rpois -> NA)
      logprof <- pmin(logprof, log(2e5))
      cnt <- rpois(n_i[d] * G, lambda = exp(logprof))
      stopifnot(!any(is.na(cnt)))
      rows[[length(rows) + 1L]] <- matrix(cnt, nrow = n_i[d], byrow = TRUE)
      meta[[length(meta) + 1L]] <- data.frame(cell_id = sprintf("%s_%s_%s", ct, d, seq_len(n_i[d])),
                                              cell_type = ct, donor = d, site = site[[d]],
                                              stringsAsFactors = FALSE)
    }
  }
  counts <- do.call(rbind, rows)
  rownames(counts) <- sprintf("c%06d", seq_len(nrow(counts)))
  colnames(counts) <- sprintf("g%04d", seq_len(G))
  meta <- do.call(rbind, meta)
  rownames(meta) <- rownames(counts)
  list(counts = counts, meta = meta, shared_genes = shared_genes, shared_z = shared_z,
       priv_genes = priv_genes, priv_z = priv_z, site = site, cts = cts)
}

## ---- anchor: Stage A v0.2.0 recipe on donor pseudobulk ---------------------
build_anchor <- function(sim, max_factors = 30, n_var = 2000, seed = 1) {
  cts <- sim$cts
  pb <- lapply(cts, function(ct) {
    sel <- sim$meta$cell_type == ct
    out <- rowsum(sim$counts[sel, , drop = FALSE], group = sim$meta$donor[sel])
    out[order(rownames(out)), , drop = FALSE]
  })
  names(pb) <- cts
  sm <- unique(sim$meta[, c("donor", "site")])
  rownames(sm) <- sm$donor
  sm$observation_id <- sm$donor
  x <- as_cell_program_data(pb, sample_metadata = sm)
  ## C4 guard: K <= n - q - 1 per cell type; clamp for pilot-scale data
  ## (the registered K30 applies at cohort scale)
  n_min <- min(vapply(pb, nrow, numeric(1)))
  q <- 1L  # one covariate column ("site")
  K <- min(max_factors, as.integer(n_min - q - 1))
  if (K < max_factors) message(sprintf("anchor K clamped to %d (pilot scale, n_min=%d)", K, n_min))
  fit <- fit_celltype_programs(x, max_factors = K,
                               loading_prior = "point_laplace", var_type = 1,
                               features = "variable", n_variable_genes = n_var,
                               covariates = "site", covariate_mode = "residualize",
                               seed = seed)
  canonicalize_programs(fit)
}

## ---- cell-level query matrix (log1p CPM per cell) ---------------------------
query_matrix <- function(sim) {
  L1 <- log1p(t(t(sim$counts) / colSums(sim$counts)) * 1e4)  # cells x genes, CPM->log1p
  L1[sim$meta$cell_type %in% sim$cts, , drop = FALSE]
}

## ---- anchored refit arms -----------------------------------------------------
## X_cell: cells x genes log-normalized query; L_anchor: genes x K anchor
## loadings (anchor$loadings[[ct]]); A0: cells x K initial activities.
## Returns list(activities, loadings, fit).
arm_anchored <- function(X_cell, L_anchor, A0, fixed_side = c("none", "activity", "effects"),
                         add_factors = 0L, seed = 1, S = NULL) {
  fixed_side <- match.arg(fixed_side)
  stopifnot(requireNamespace("flashier.utils", quietly = TRUE))
  f <- flashier.utils::fit_initialized_ebmf(X = X_cell, sample_side = "rows",
                            activity_init = A0, effect_init = L_anchor,
                            activity_prior = ebnm::ebnm_point_normal,
                            effect_prior = ebnm::ebnm_point_laplace,
                            fixed_side = fixed_side, add_factors = as.integer(add_factors),
                            var_type = 1L, feature_scale = "logcounts",
                            S = S, seed = seed, verbose = 0L)
  .fu_extract(f, sample_dim = "rows")
}

## flashier.utils fit -> list(activities = samples x factors, effects = features x
## factors), via standardize_factors (the package's own extraction convention,
## R/ebmf.R); falls back to the raw fields when the metadata attr is absent.
.fu_extract <- function(f, sample_dim = "rows") {
  if (!is.null(attr(f, "flashbridge_metadata"))) {
    rms <- flashier.utils::standardize_factors(f, sample_side = sample_dim,
                                               scaling = "rms", orientation = "anchor_feature")
    list(activities = as.matrix(rms$activity), effects = as.matrix(rms$effects), fit = f)
  } else {
    list(activities = as.matrix(f$activity %||% f$Z), effects = as.matrix(f$effects %||% f$W), fit = f)
  }
}

## initialize cell activities from anchor loadings: plain OLS on the query
## (the real-data variant uses project_programs_to_cells() instead).
ols_activity <- function(X_cell, L) {
  A <- t(solve(crossprod(L), crossprod(L, t(X_cell))))   # cells x K
  dimnames(A) <- list(rownames(X_cell), colnames(L))
  A
}

## named donor lookup: meta$donor indexed by cell id (an UNNAMED vector
## indexed by character names returns all-NA)
donor_of <- function(meta) setNames(meta$donor, rownames(meta))

## metacell grouting for the mc arm: k-means within donor, sum counts,
## S = sqrt(size) (Poisson-scale SEs on the summed counts).
metacell_grout <- function(counts_ct, meta_ct, n_mc_per_donor = 20, seed = 1) {
  set.seed(seed)
  Lg <- log1p(t(t(counts_ct) / colSums(counts_ct)) * 1e4)   # cells x genes
  donors <- donor_of(meta_ct)[rownames(Lg)]
  parts <- lapply(split(seq_len(nrow(Lg)), donors), function(idx) {
    k <- min(n_mc_per_donor, length(idx))
    km <- kmeans(Lg[idx, , drop = FALSE], centers = k, iter.max = 30, nstart = 1)
    data.frame(cell = rownames(Lg)[idx], mc = km$cluster, donor = donors[idx][1],
               stringsAsFactors = FALSE)
  })
  out <- do.call(rbind, parts)
  out$mc_id <- sprintf("%s_mc%02d", out$donor, out$mc)
  out
}

## ---- fit-contract wrapper for the CMP caller --------------------------------
## find_coordinated_programs() consumes fit$scores / fit$loadings /
## fit$observation_ids / fit$sample_metadata / fit$cell_types.
wrap_for_caller <- function(scores, loadings, donors, sample_metadata, cts) {
  structure(list(scores = scores, loadings = loadings, observation_ids = donors,
                 sample_metadata = sample_metadata, cell_types = cts),
            class = c("cell_program_fit", "list"))
}

## ---- metrics -----------------------------------------------------------------
## drift: Hungarian-matched mean |cor(anchor donor scores, refit donor-mean
## activities)| over anchor programs; reported as 1 - value.
drift_vs_anchor <- function(anchor_scores, refit_scores) {
  common <- intersect(rownames(anchor_scores), rownames(refit_scores))
  if (length(common) < 3L) return(NA_real_)
  A <- scale(anchor_scores[common, , drop = FALSE])
  B <- scale(refit_scores[common, , drop = FALSE])
  C <- abs(stats::cor(A, B)); C[is.na(C)] <- 0
  if (ncol(C) < 1L) return(NA_real_)
  k <- min(ncol(A), ncol(B))
  sol <- clue::solve_LSAP(C, maximum = TRUE)
  mean(C[cbind(seq_len(k), sol[seq_len(k)])])
}

## per-cell activity stability: n_sub subsamples of `frac` of cells; re-project
## the arm's loadings onto the subsample; mean Spearman on shared cells.
activity_stability <- function(X_cell, activities, loadings, n_sub = 10L, frac = 0.8, seed = 1) {
  set.seed(seed)
  v <- apply(activities, 2, stats::var)
  keep_k <- names(v)[v > 1e-12]
  if (!length(keep_k)) return(NA_real_)
  L <- as.matrix(loadings[, keep_k, drop = FALSE])
  A0 <- activities[, keep_k, drop = FALSE]
  vals <- replicate(n_sub, {
    idx <- sample(nrow(X_cell), max(10L, floor(frac * nrow(X_cell))))
    A2 <- ols_activity(X_cell[idx, , drop = FALSE], L)
    mean(vapply(seq_along(keep_k), function(j)
      suppressWarnings(stats::cor(A0[rownames(A2), j], A2[, j])), numeric(1)), na.rm = TRUE)
  })
  mean(unlist(vals), na.rm = TRUE)
}

## reconstruction guardrail: mean per-gene residual variance ratio vs denovo
## (ELBO-equivalent proxy; > 1.01 flags the guardrail).
recon_ratio <- function(X_cell, activities, loadings, ref_activities, ref_loadings) {
  .mrv <- function(A, L) {
    R <- X_cell - A %*% t(L)
    mean(apply(R, 1, stats::var, na.rm = TRUE))
  }
  .mrv(activities, loadings) / .mrv(ref_activities, ref_loadings)
}
