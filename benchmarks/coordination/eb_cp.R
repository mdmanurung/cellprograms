## eb_cp.R — empirical-Bayes CP tensor factorization (genes x donors x cell types)
##
## Model: Y_ct[d, g] = sum_r U[g, r] V[d, r] W[ct, r] + eps, over observed
## (ct, d) pairs only (ragged donor x cell-type grid is native; unobserved
## slices contribute nothing). Flash-style fitting: greedy rank-1 components,
## each mode update is a normal-means problem solved by EBNM with
## per-component estimated priors (point_laplace on genes and cell types,
## normal on donors), then backfitting.
##
## Variance model: per-(cell-type, donor) residual variance tau2[[ct]][d]
## (the tensor analogue of flashier var_type = 1 per-row variance, the
## package's validated pseudobulk default — pseudobulk noise is heteroskedastic
## across donors via cell counts, and a per-cell-type variance mis-calibrates
## every normal-means SE and the held-out statistic).
##
## Model selection (FROZEN): winner's-curse-proof held-out component test.
## Donors are split 50/50 synchronously across cell types, stratified by
## `strata`. Components are greedily fit on half A; component r is kept iff
## its held-out SS reduction on half B (free per-row donor coefficients,
## fixed A-side gene/cell-type modes, per-row studentization) exceeds the
## noise scale:
##   z_r = (sum_{ct,d in B} (Y_ct[d,].u)^2 / (tau2_ct[d] ||u||^2) - n_B) / sqrt(2 n_B)
## kept iff z_r > heldout_z (= 3), stopping at the first failure. Kept
## components are refit on the full data, warm-started, and backfit.
##
## Identifiability: CP components are essentially unique (Kruskal), so the
## cell-type mode W needs no post-hoc matching. Membership rule (FROZEN):
## ct is a member of component r iff |W[ct, r]| >= member_thresh * max_c
## |W[c, r]|, member_thresh = 0.2. A component with >= 2 members is a
## coordination call.

## SVD seed for a rank-1 step (flash-style greedy init): strongest per-cell-type
## top singular triplet (row-standardized), embedded in the union donor space.
.cp_seed <- function(R, obs, tau2) {
  best <- list(s = -Inf)
  for (ci in seq_along(R)) {
    Rs <- R[[ci]] / sqrt(tau2[[ci]])              # whiten rows for the seed
    sv <- svd(Rs, nu = 1L, nv = 1L)
    if (sv$d[1] > best$s) best <- list(s = sv$d[1], ci = ci, u = sv$v[, 1], v = sv$u[, 1])
  }
  v <- numeric(length(obs$donors)); v[obs$idx[[best$ci]]] <- best$v
  w <- numeric(length(R)); w[best$ci] <- 1
  list(u = best$u, v = v, w = w)
}

.cp_rank1 <- function(R, obs, tau2, ebnm_gene, ebnm_donor, ebnm_ct,
                      init = NULL, max_iter = 100L, tol = 1e-6) {
  C <- length(R); G <- ncol(R[[1]]); D <- length(obs$donors)
  if (is.null(init)) init <- .cp_seed(R, obs, tau2)
  v <- init$v; w <- init$w; u <- init$u
  prev_obj <- -Inf
  it <- 0L
  for (it in seq_len(max_iter)) {
    ## ---- gene mode ----
    num <- numeric(G); den <- 0; s2 <- 0
    for (ci in seq_len(C)) {
      ii <- obs$idx[[ci]]; vv <- v[ii]
      num <- num + w[ci] * drop(crossprod(vv, R[[ci]]))
      s2v_t <- sum(vv^2 * tau2[[ci]])            # sum_d v_d^2 tau2_ct[d]
      den <- den + w[ci]^2 * sum(vv^2)
      s2 <- s2 + w[ci]^2 * s2v_t
    }
    if (den < 1e-12) return(NULL)
    u <- as.numeric(ebnm_gene(num / den, sqrt(s2) / den)$posterior$mean)
    u2 <- sum(u^2)
    if (u2 < 1e-12) return(NULL)
    ## ---- donor mode ----
    num <- numeric(D); den <- numeric(D); s2 <- numeric(D)
    for (ci in seq_len(C)) {
      ii <- obs$idx[[ci]]
      proj <- drop(R[[ci]] %*% u)
      num[ii] <- num[ii] + w[ci] * proj
      den[ii] <- den[ii] + w[ci]^2 * u2
      s2[ii] <- s2[ii] + w[ci]^2 * u2 * tau2[[ci]]
    }
    ok <- den > 1e-12
    if (!any(ok)) return(NULL)
    v[ok] <- as.numeric(ebnm_donor(num[ok] / den[ok], sqrt(s2[ok]) / den[ok])$posterior$mean)
    v[!ok] <- 0
    ## ---- cell-type mode ----
    num <- numeric(C); den <- numeric(C); s2 <- numeric(C)
    for (ci in seq_len(C)) {
      ii <- obs$idx[[ci]]; vv <- v[ii]
      num[ci] <- sum(vv * drop(R[[ci]] %*% u))
      den[ci] <- sum(vv^2) * u2
      s2[ci] <- sum(vv^2 * tau2[[ci]]) * u2
    }
    okc <- den > 1e-12
    if (!any(okc)) return(NULL)
    w[okc] <- as.numeric(ebnm_ct(num[okc] / den[okc], sqrt(s2[okc]) / den[okc])$posterior$mean)
    w[!okc] <- 0
    ## ---- per-(ct, donor) residual variance + objective ----
    ss <- 0
    for (ci in seq_len(C)) {
      E <- R[[ci]] - w[ci] * tcrossprod(v[obs$idx[[ci]]], u)
      tau2[[ci]] <- pmax(rowMeans(E^2), 1e-8)
      ss <- ss + sum(E^2)
    }
    obj <- -ss / 2
    if (abs(obj - prev_obj) < tol * (1 + abs(obj))) break
    prev_obj <- obj
  }
  if (sqrt(sum(u^2)) < 1e-6 || max(abs(w)) < 1e-6) return(NULL)
  list(u = u, v = v, w = w, tau2 = tau2, n_iter = it)
}

## held-out statistic for one component: SS reduction on B with fixed (u, w),
## free per-row donor coefficients, per-row studentization.
.cp_heldout_z <- function(u, w, mats_B, tau2_B) {
  T_stat <- 0; n_B <- 0
  u2 <- sum(u^2)
  if (u2 < 1e-12) return(-Inf)
  for (ci in seq_along(mats_B)) {
    Y <- mats_B[[ci]]
    if (nrow(Y) == 0L) next
    proj <- drop(Y %*% u)
    T_stat <- T_stat + sum(proj^2 / (tau2_B[[ci]] * u2))
    n_B <- n_B + nrow(Y)
  }
  (T_stat - n_B) / sqrt(2 * n_B)
}

.split_donors <- function(donors, strata, seed) {
  set.seed(seed)
  inA <- rep(FALSE, length(donors)); names(inA) <- donors
  st <- if (is.null(strata)) rep("all", length(donors)) else strata[donors]
  for (g in unique(st)) {
    d <- donors[st == g]
    inA[sample(d, floor(length(d) / 2))] <- TRUE
  }
  inA
}

## per-row variance of a matrix, floored
.rowvar <- function(m) pmax(matrixStats::rowVars(m), 1e-8)

#' EB-CP tensor factorization
#'
#' @param mats Named list of (donors_ct x genes) matrices, common gene columns,
#'   rownames = donor ids. Already covariate-adjusted and gene-centered.
#' @param Kmax Maximum greedy components.
#' @param ebnm_gene,ebnm_donor,ebnm_ct EBNM functions per mode.
#' @param member_thresh Frozen membership threshold on max-normalized |W|.
#' @param heldout_z Held-out z threshold for keeping a component (frozen: 3).
#' @param strata Named donor strata for the split (e.g. institute).
#' @param backfit_iter Backfitting sweeps after the greedy pass.
#' @param seed Random seed.
eb_cp_fit <- function(mats, Kmax = 15L,
                      ebnm_gene = ebnm::ebnm_point_laplace,
                      ebnm_donor = ebnm::ebnm_normal,
                      ebnm_ct = ebnm::ebnm_point_laplace,
                      member_thresh = 0.2, heldout_z = 3,
                      strata = NULL, max_iter = 100L, tol = 1e-6,
                      backfit_iter = 2L, seed = 1L) {
  cts <- names(mats)
  genes <- Reduce(intersect, lapply(mats, colnames))
  mats <- lapply(mats, function(m) m[, genes, drop = FALSE])
  donors <- sort(unique(unlist(lapply(mats, rownames))))
  obs <- list(donors = donors,
              idx = lapply(mats, function(m) match(rownames(m), donors)))

  ## ---- donor-synchronous split, stratified ----
  inA <- .split_donors(donors, strata, seed)
  mats_A <- lapply(mats, function(m) m[rownames(m) %in% donors[inA], , drop = FALSE])
  mats_B <- lapply(mats, function(m) m[rownames(m) %in% donors[!inA], , drop = FALSE])
  obs_A <- list(donors = donors[inA],
                idx = lapply(mats_A, function(m) match(rownames(m), donors[inA])))
  tau2_A <- lapply(mats_A, .rowvar)

  ## ---- greedy on half A with held-out selection on half B ----
  R_A <- lapply(mats_A, function(m) m + 0)
  comps_A <- list(); kept <- numeric(0); trace_z <- numeric(0)
  set.seed(seed)
  for (r in seq_len(Kmax)) {
    f <- .cp_rank1(R_A, obs_A, tau2_A, ebnm_gene, ebnm_donor, ebnm_ct,
                   max_iter = max_iter, tol = tol)
    if (is.null(f)) break
    ## B-side per-row studentization: each held-out row's own marginal variance
    ## (per-donor noise is not predictable for unseen donors; the row's marginal
    ## variance estimates it directly, at negligible O(1/G) self-dependence).
    tau2_B <- lapply(mats_B, function(mb) if (nrow(mb)) .rowvar(mb) else numeric(0))
    z <- .cp_heldout_z(f$u, f$w, mats_B, tau2_B)
    trace_z <- c(trace_z, z)
    if (!is.finite(z) || z <= heldout_z) break
    comps_A[[length(comps_A) + 1L]] <- f
    kept <- c(kept, z)
    for (ci in seq_along(cts)) {
      R_A[[ci]] <- R_A[[ci]] - f$w[ci] * tcrossprod(f$v[obs_A$idx[[ci]]], f$u)
    }
    tau2_A <- f$tau2
  }

  ## ---- refit kept components on full data, warm-started, then backfit ----
  set.seed(seed + 1L)
  tau2 <- lapply(mats, .rowvar)
  R <- lapply(mats, function(m) m + 0)
  comps <- list()
  for (r in seq_along(comps_A)) {
    init <- comps_A[[r]]
    v_full <- numeric(length(donors))                 # embed A-side donor mode
    v_full[match(donors[inA], donors)] <- init$v      # B donors start at 0
    init$v <- v_full
    f <- .cp_rank1(R, obs, tau2, ebnm_gene, ebnm_donor, ebnm_ct,
                   init = init, max_iter = max_iter, tol = tol)
    if (is.null(f)) next
    comps[[length(comps) + 1L]] <- f
    for (ci in seq_along(cts)) {
      R[[ci]] <- R[[ci]] - f$w[ci] * tcrossprod(f$v[obs$idx[[ci]]], f$u)
    }
    tau2 <- f$tau2
  }
  if (length(comps) > 1L) {
    for (sweep in seq_len(backfit_iter)) {
      drop <- integer(0)
      for (r in seq_along(comps)) {
        f <- comps[[r]]
        for (ci in seq_along(cts)) {
          R[[ci]] <- R[[ci]] + f$w[ci] * tcrossprod(f$v[obs$idx[[ci]]], f$u)
        }
        f2 <- .cp_rank1(R, obs, tau2, ebnm_gene, ebnm_donor, ebnm_ct,
                        init = f, max_iter = max_iter, tol = tol)
        if (is.null(f2)) { drop <- c(drop, r); next }
        comps[[r]] <- f2
        for (ci in seq_along(cts)) {
          R[[ci]] <- R[[ci]] - f2$w[ci] * tcrossprod(f2$v[obs$idx[[ci]]], f2$u)
        }
        tau2 <- f2$tau2
      }
      if (length(drop)) comps <- comps[-drop]
    }
  }
  R_tot <- sum(vapply(R, function(m) sum(m^2), numeric(1)))
  Y_tot <- sum(vapply(mats, function(m) sum(m^2), numeric(1)))

  R_n <- length(comps)
  U <- matrix(0, length(genes), R_n, dimnames = list(genes, sprintf("CP%02d", seq_len(R_n))))
  V <- matrix(0, length(donors), R_n, dimnames = list(donors, colnames(U)))
  W <- matrix(0, length(cts), R_n, dimnames = list(cts, colnames(U)))
  for (r in seq_len(R_n)) {
    f <- comps[[r]]
    sgn <- sign(f$u[which.max(abs(f$u))]); if (sgn == 0) sgn <- 1
    sc <- sqrt(mean(f$v^2)); if (sc < 1e-12) sc <- 1
    U[, r] <- f$u * sgn * sc
    V[, r] <- f$v * sgn / sc
    W[, r] <- f$w * sgn
  }
  list(U = U, V = V, W = W, genes = genes, donors = donors, cts = cts,
       tau2 = tau2, member_thresh = member_thresh, heldout_z = heldout_z,
       kept_heldout_z = kept, trace_z = trace_z,
       frac_unexplained = R_tot / Y_tot, n_components = R_n)
}

#' Coordination calls from an eb_cp fit (benchmark output contract)
eb_cp_clusters <- function(fit, member_thresh = fit$member_thresh) {
  empty <- data.frame(cluster = integer(), cell_type = character(),
                      program_id = character(), stat = numeric(),
                      fwer_p = numeric(), stringsAsFactors = FALSE)
  if (fit$n_components == 0L) return(list(clusters = empty, activity = NULL))
  rows <- list(); acts <- list(); id <- 0L
  for (r in seq_len(fit$n_components)) {
    w <- fit$W[, r]
    mem <- names(w)[abs(w) >= member_thresh * max(abs(w))]
    if (length(mem) < 2L) next
    id <- id + 1L
    rows[[id]] <- data.frame(cluster = id, cell_type = mem,
                             program_id = colnames(fit$W)[r],
                             stat = max(abs(w)), fwer_p = NA_real_,
                             stringsAsFactors = FALSE)
    acts[[id]] <- fit$V[, r]
  }
  cl <- if (length(rows)) do.call(rbind, rows) else empty
  activity <- NULL
  if (length(acts)) {
    activity <- do.call(cbind, acts)
    colnames(activity) <- paste0("cluster_", seq_along(acts))
    rownames(activity) <- fit$donors
  }
  list(clusters = cl, activity = activity)
}
