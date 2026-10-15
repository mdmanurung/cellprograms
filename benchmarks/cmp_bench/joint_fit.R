## joint_fit.R -- PLAN v2 option 2: MOFA+ joint factorization across cell types, cross-fitted per-cell-type scores.
## Rscript joint_fit.R <base> <scenarios,comma> <seed_from> <seed_to> <amp> <cores> [K=15] [outroot]
## Writes two cache objects per dataset (same layout as fit_cache.R):
##   fits/jointA/<base>/<sc>/sNNNN.rds  : joint (factor, view) nodes only (cross-fitted hybrid, J-A)
##   fits/jointB/<base>/<sc>/sNNNN.rds  : main-arm per-cell-type programs + joint nodes (augment, J-B)
## Input per view = same preprocessing as fit_celltype_programs (2000 HVG, residualize on institute, center).
## Cross-fitting: 5 site-stratified donor folds; MOFA fit on 4/5; held-out donors scored by per-view least-squares
## projection onto W_v; fold factors matched to the full-data fit by |cor| of factor scores on training donors (Hungarian), a factor
## is kept only if matched (|r| >= 0.5) in all folds. A node (k, v) is stored with its full-fit view-v R2 (node_r2); R2 filtering happens in the candidate (cand_v2.R).
source("/mnt/shared-workspace/cmp/bench/common.R")
suppressPackageStartupMessages(library(parallel))

prep_views <- function(x, n_var = 2000L, covariates = "institute") {
  fs <- .select_features(x, "variable", n_var)
  V <- lapply(stats::setNames(x$cell_types, x$cell_types), function(ct) {
    Y <- as.matrix(x$pseudobulk[[ct]])[, fs[[ct]], drop = FALSE]
    Y <- .residualize(Y, x$sample_metadata, covariates, ct)
    Y <- scale(Y, center = TRUE, scale = FALSE)
    Y[, colSums(abs(Y)) > 1e-8 * max(colSums(abs(Y)), 1), drop = FALSE]
  })
  V[vapply(V, function(Y) nrow(Y) >= 20 && ncol(Y) >= 50, TRUE)]
}

stack_views <- function(V, donors) lapply(V, function(Y) { M <- matrix(NaN, length(donors), ncol(Y), dimnames = list(donors, colnames(Y)))
  M[intersect(donors, rownames(Y)), ] <- Y[intersect(donors, rownames(Y)), ]; M })

run_mofa <- function(V, donors, K, seed) {
  mh <- reticulate::import_from_path("mofa_helper", "/mnt/shared-workspace/cmp/bench")
  r <- mh$fit_mofa(unname(lapply(stack_views(V, donors), function(M) M)), K = as.integer(K), seed = as.integer(seed))
  W <- lapply(seq_along(V), function(m) { w <- as.matrix(r$W[[m]]); rownames(w) <- colnames(V[[m]]); w }); names(W) <- names(V)
  list(W = W, Z = as.matrix(r$Z), r2 = as.matrix(r$r2))
}

proj_scores <- function(Y, W) {   # least-squares projection of rows of Y onto columns of W
  g <- intersect(colnames(Y), rownames(W)); Wg <- W[g, , drop = FALSE]
  A <- crossprod(Wg); A <- A + diag(1e-6 * max(sum(diag(A)), 1e-12), ncol(A))
  Y[, g, drop = FALSE] %*% Wg %*% solve(A)
}

joint_one <- function(base, sc, seed, amp, K, outroot) {
  fA <- file.path(outroot, "jointA", base, sc, sprintf("s%04d.rds", seed)); fB <- sub("/jointA/", "/jointB/", fA)
  if (file.exists(fA) && file.size(fA) > 0 && file.exists(fB) && file.size(fB) > 0) return("skip")
  t0 <- Sys.time()
  d <- gen_cmp(base, sc, seed, amp)
  V <- prep_views(d$x)
  donors <- sort(unique(unlist(lapply(V, rownames))))
  full <- run_mofa(V, donors, K, seed)
  Kf <- ncol(full$Z)
  site <- d$x$sample_metadata[donors, "institute"]
  set.seed(seed + 99L)
  fold <- stats::setNames(integer(length(donors)), donors)
  for (s in unique(site)) { w <- which(site == s); fold[w] <- sample(rep_len(1:5, length(w))) }
  stackW <- function(W) do.call(rbind, lapply(names(W), function(v) { w <- W[[v]]; rownames(w) <- paste(v, rownames(w)); w }))
  SW <- stackW(full$W)
  S <- lapply(V, function(Y) matrix(NA_real_, nrow(Y), Kf, dimnames = list(rownames(Y), NULL)))
  matched <- rep(TRUE, Kf)
  for (f in 1:5) {
    tr <- donors[fold != f]; ho <- donors[fold == f]
    Vt <- lapply(V, function(Y) Y[intersect(rownames(Y), tr), , drop = FALSE])
    ff <- run_mofa(Vt, tr, K, seed + f)
    ## match fold factors to full-fit factors by factor scores on the training donors (loadings are too noisy)
    C <- suppressWarnings(stats::cor(full$Z[match(tr, donors), , drop = FALSE], ff$Z)); C[!is.finite(C)] <- 0
    nf <- ncol(C); pad <- max(Kf, nf); Cp <- matrix(0, pad, pad); Cp[seq_len(Kf), seq_len(nf)] <- abs(C)
    asg <- as.integer(clue::solve_LSAP(Cp, maximum = TRUE))[seq_len(Kf)]
    for (k in seq_len(Kf)) {
      j <- asg[k]
      if (j > nf || abs(C[k, j]) < 0.5) { matched[k] <- FALSE; next }
      sg <- sign(C[k, j])
      for (v in names(V)) {
        h <- intersect(rownames(V[[v]]), ho); if (!length(h)) next
        S[[v]][h, k] <- sg * proj_scores(V[[v]][h, , drop = FALSE], ff$W[[v]])[, j]
      }
    }
  }
  keepK <- which(matched)
  scores <- list(); loadings <- list(); node_r2 <- c()
  for (m in seq_along(V)) {
    v <- names(V)[m]
    ks <- keepK
    ks <- ks[colSums(is.na(S[[v]][, ks, drop = FALSE])) == 0]
    if (!length(ks)) next
    nm <- paste0(v, "_J", ks)
    scores[[v]] <- S[[v]][, ks, drop = FALSE]; colnames(scores[[v]]) <- nm
    loadings[[v]] <- full$W[[v]][, ks, drop = FALSE]; colnames(loadings[[v]]) <- nm
    node_r2[nm] <- full$r2[ks, m]
  }
  rt <- as.numeric(Sys.time() - t0, units = "secs")
  fitA <- list(cell_types = names(scores), observation_ids = donors, scores = scores, loadings = loadings,
               joint = list(K = Kf, kept = keepK, r2 = full$r2, view_names = names(V), node_r2 = node_r2))
  objA <- list(fit = fitA, truth = d$truth, strata = d$strata, context = d$context, ncell = d$ncell, meta = d$x$sample_metadata,
               arm = "jointA", amp = amp, runtime = rt)
  ## J-B: augment main-arm programs with the joint nodes
  fm <- file.path(CMP_ROOT, "fits/main", base, sc, sprintf("s%04d.rds", seed))
  objB <- NULL
  if (file.exists(fm)) {
    om <- readRDS(fm)
    if (!inherits(om$fit, "fit_error")) {
      fb <- om$fit
      for (v in names(scores)) {
        if (is.null(fb$scores[[v]]) || !ncol(fb$scores[[v]])) { fb$scores[[v]] <- scores[[v]]; fb$loadings[[v]] <- loadings[[v]]; next }
        ids <- intersect(rownames(fb$scores[[v]]), rownames(scores[[v]]))
        g <- intersect(rownames(fb$loadings[[v]]), rownames(loadings[[v]]))
        fb$scores[[v]] <- cbind(fb$scores[[v]][ids, , drop = FALSE], scores[[v]][ids, , drop = FALSE])
        fb$loadings[[v]] <- cbind(fb$loadings[[v]][g, , drop = FALSE], loadings[[v]][g, , drop = FALSE])
      }
      fb$cell_types <- union(fb$cell_types, names(scores))
      fb$joint <- fitA$joint; objB <- modifyList(om, list(fit = fb, arm = "jointB", runtime = rt))
    }
  }
  sv <- function(obj, f) { if (is.null(obj)) return(invisible()); dir.create(dirname(f), recursive = TRUE, showWarnings = FALSE)
    tmp <- tempfile(fileext = ".rds"); saveRDS(obj, tmp); ok <- system2("cp", c(shQuote(tmp), shQuote(f))) == 0L; unlink(tmp)
    if (!ok || file.size(f) == 0) stop("copy failed: ", f) }
  sv(objA, fA); sv(objB, fB)
  sprintf("%.0fs K=%d kept=%d nodes=%d", rt, Kf, length(keepK), sum(vapply(scores, ncol, 1L)))
}

if (sys.nframe() == 0L) {
  a <- commandArgs(TRUE)
  base <- a[1]; scs <- strsplit(a[2], ",")[[1]]; s0 <- as.integer(a[3]); s1 <- as.integer(a[4]); amp <- as.numeric(a[5])
  cores <- as.integer(a[6]); K <- if (length(a) > 6) as.integer(a[7]) else 15L
  outroot <- if (length(a) > 7) a[8] else file.path(CMP_ROOT, "fits")
  jobs <- expand.grid(seed = s0:s1, sc = scs, stringsAsFactors = FALSE)
  res <- mclapply(seq_len(nrow(jobs)), function(i) {
    r <- tryCatch(joint_one(base, jobs$sc[i], jobs$seed[i], amp, K, outroot), error = function(e) paste("ERR", conditionMessage(e)))
    cat(sprintf("joint %s %s s%d %s\n", base, jobs$sc[i], jobs$seed[i], r)); r
  }, mc.cores = cores, mc.preschedule = FALSE)
  cat("finished", sum(grepl("ERR", unlist(res))), "errors\n")
}
