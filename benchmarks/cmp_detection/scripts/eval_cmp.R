## eval_cmp.R -- frozen evaluator for coordinated multicellular programs (cmp-bench-v1)
## Node = (cell_type, program). A node belongs to the signal (target/decoy/context/SITE) with the
## largest |cor(score, z_ct)| among signals that include its cell type, if that |cor| >= NODE_THR.
## A called CMP is labelled by the majority signal of its nodes ("spurious" if majority unassigned).
## F1(signal) = max over CMPs of harmonic mean(recall = member cts with an assigned node in the CMP / |members|,
##                                             precision = nodes assigned to signal / CMP size).
NODE_THR <- 0.4
POS_SCEN <- c("cmp_distinct", "cmp_mixed", "context", "ambient", "lr", "weak")

## within-context permutation strata: site x context level (context binned to quintiles if > 5 levels).
## Residualizing on context removes shared means but not shared context-dependent variance, which makes
## donors non-exchangeable across context levels; permuting within context restores exchangeability.
within_strata <- function(strata, context) {
  if (is.null(context)) return(strata)
  cx <- context[names(strata)]
  lev <- if (length(unique(cx)) > 5) cut(cx, unique(stats::quantile(cx, 0:5 / 5)), include.lowest = TRUE) else cx
  stats::setNames(paste(strata, lev, sep = ":"), names(strata))
}
## within-context transform: residualize each cell type's scores on context (+ intercept)
within_context <- function(fit, context) {
  if (is.null(context)) return(fit)
  fit$scores <- lapply(fit$scores, function(S) {
    if (!ncol(S)) return(S)
    x <- context[rownames(S)]; X <- cbind(1, x)
    R <- S - X %*% solve(crossprod(X), crossprod(X, S)); dimnames(R) <- dimnames(S); R
  })
  fit
}

.signal_table <- function(truth, strata, mode) {
  sg <- truth$signals
  typ <- vapply(sg, function(s) if (s$type == "context") (if (mode == "within") "decoy" else "target") else s$type, "")
  ## SITE pseudo-signal: standardized site indicator projections (all cell types); a node is SITE if R2(site) >= thr^2
  list(sg = sg, type = setNames(typ, names(sg)))
}

node_assign <- function(fit, truth, strata, mode = "marginal") {
  st <- .signal_table(truth, strata, mode)
  out <- list()
  for (ct in fit$cell_types) {
    S <- fit$scores[[ct]]; if (is.null(S) || !ncol(S)) next
    d <- rownames(S)
    cand <- names(st$sg)[vapply(st$sg, function(s) ct %in% s$members, TRUE)]
    C <- matrix(0, ncol(S), length(cand) + 1L, dimnames = list(colnames(S), c(cand, "SITE")))
    for (s in cand) { z <- st$sg[[s]]$z[[ct]][d]; C[, s] <- abs(suppressWarnings(stats::cor(S, z, use = "pairwise"))) }
    f <- factor(strata[d])
    C[, "SITE"] <- if (nlevels(f) > 1) sqrt(pmax(0, apply(S, 2, function(v) summary(stats::lm(v ~ f))$r.squared))) else 0
    C[is.na(C)] <- 0
    best <- colnames(C)[max.col(C, ties.method = "first")]; bv <- apply(C, 1, max)
    out[[ct]] <- data.frame(cell_type = ct, program_id = colnames(S), signal = ifelse(bv >= NODE_THR, best, NA),
                            cor = bv, stringsAsFactors = FALSE)
  }
  na <- do.call(rbind, out)
  na$type <- ifelse(is.na(na$signal), NA, ifelse(na$signal == "SITE", "site", st$type[na$signal]))
  na
}

## evaluate one called-cluster table against truth
eval_calls <- function(cl, fit, truth, strata, mode = "marginal") {
  na <- node_assign(fit, truth, strata, mode)
  st <- .signal_table(truth, strata, mode)
  tgt <- names(st$type)[st$type == "target"]
  ct_level <- nrow(cl) > 0 && all(is.na(cl$program_id))
  key <- paste(na$cell_type, na$program_id)
  cmps <- if (nrow(cl)) split(cl, cl$cluster) else list()
  lab <- vapply(cmps, function(c) {
    if (ct_level) {   # cell-type-level clusters: label by signal covering most cts (a ct counts if any node assigned)
      sigs <- unique(na$signal[!is.na(na$signal) & na$cell_type %in% c$cell_type])
      if (!length(sigs)) return("spurious")
      cov <- vapply(sigs, function(s) length(intersect(c$cell_type, na$cell_type[na$signal %in% s])), 1)
      return(if (max(cov) / nrow(c) > 0.5) sigs[which.max(cov)] else "spurious")
    }
    s <- na$signal[match(paste(c$cell_type, c$program_id), key)]
    tb <- table(s, useNA = "no")
    if (!length(tb) || max(tb) / nrow(c) <= 0.5) "spurious" else names(tb)[which.max(tb)]
  }, "")
  lab_type <- ifelse(lab == "spurious", "spurious", ifelse(lab == "SITE", "site", st$type[lab]))
  f1 <- vapply(tgt, function(s) {
    mem <- st$sg[[s]]$members
    if (!length(cmps)) return(0)
    max(vapply(cmps, function(c) {
      if (ct_level) {
        hit <- intersect(c$cell_type, intersect(mem, na$cell_type[na$signal %in% s]))
        rec <- length(hit) / length(mem); pre <- length(intersect(c$cell_type, mem)) / nrow(c)
      } else {
        s_n <- na$signal[match(paste(c$cell_type, c$program_id), key)]
        ok <- !is.na(s_n) & s_n == s
        rec <- length(unique(c$cell_type[ok])) / length(mem); pre <- mean(ok)
      }
      if (rec + pre == 0) 0 else 2 * rec * pre / (rec + pre)
    }, 1))
  }, 1)
  oracle <- vapply(tgt, function(s) length(unique(na$cell_type[na$signal %in% s])) / length(st$sg[[s]]$members), 1)
  list(f1 = f1, oracle = oracle, labels = lab, label_type = lab_type,
       summary = data.frame(n_cmp = length(cmps), n_false = sum(lab_type != "target"),
                            any_false = any(lab_type != "target"),
                            any_decoy = any(lab_type == "decoy"), any_site = any(lab_type == "site"),
                            any_spurious = any(lab_type == "spurious"),
                            mean_f1 = if (length(f1)) mean(f1) else NA_real_,
                            mean_oracle = if (length(oracle)) mean(oracle) else NA_real_,
                            n_nodes_called = nrow(cl), stringsAsFactors = FALSE))
}

## run a method function over cached fits. method_fun(fit, strata, seed) -> clusters df
## returns per (dataset, mode) summary rows and per-signal F1 rows
eval_method <- function(method_fun, files, cores = 16L, label = "m", seed = 1L) {
  res <- parallel::mclapply(files, function(f) {
    o <- readRDS(f)
    if (inherits(o$fit, "fit_error")) return(NULL)
    modes <- if (o$truth$scenario == "context") c("marginal", "within") else "marginal"
    lapply(modes, function(md) {
      fit <- if (md == "within") within_context(o$fit, o$context) else o$fit
      if (md == "within") { fit$scores_raw <- o$fit$scores; fit$context <- o$context }
      if (exists("get_profiles")) fit$profiles <- get_profiles(o)   # method input: per-ct mean profiles (F8)
      t0 <- Sys.time()
      st <- if (md == "within") within_strata(o$strata, o$context) else o$strata
      cl <- tryCatch(method_fun(fit, st, seed), error = function(e) structure(conditionMessage(e), class = "err"))
      rt <- as.numeric(Sys.time() - t0, units = "secs")
      if (!inherits(cl, "err") && !is.data.frame(cl)) cl <- cl$clusters
      if (inherits(cl, "err")) return(list(summary = data.frame(file = f, base = o$truth$base, scenario = o$truth$scenario,
                                                                 seed = o$truth$seed, mode = md, error = as.character(cl), runtime = rt)))
      ev <- eval_calls(cl, fit, o$truth, o$strata, md)
      list(summary = cbind(data.frame(file = f, base = o$truth$base, scenario = o$truth$scenario, seed = o$truth$seed,
                                      mode = md, error = NA_character_, runtime = rt), ev$summary),
           signals = if (length(ev$f1)) data.frame(base = o$truth$base, scenario = o$truth$scenario, seed = o$truth$seed,
                                                   mode = md, signal = names(ev$f1), f1 = ev$f1, oracle = ev$oracle) else NULL,
           calls = if (nrow(cl)) cbind(data.frame(base = o$truth$base, scenario = o$truth$scenario, seed = o$truth$seed, mode = md),
                                       cl, label = unname(ev$labels[as.character(cl$cluster)])) else NULL)
    })
  }, mc.cores = cores, mc.preschedule = FALSE)
  flat <- unlist(res, recursive = FALSE)
  bind <- function(k) data.table::rbindlist(lapply(flat, `[[`, k), fill = TRUE)
  list(summary = bind("summary"), signals = bind("signals"), calls = bind("calls"), method = label)
}

## scorecard: primary metric + guardrails
scorecard <- function(ev) {
  s <- as.data.frame(ev$summary); g <- as.data.frame(ev$signals)
  if (!"any_false" %in% names(s)) s$any_false <- NA
  prim <- g[g$scenario %in% POS_SCEN & g$mode == "marginal", ]
  cp_up <- function(x) if (!length(x)) NA else stats::binom.test(sum(x), length(x))$conf.int[2]
  out <- list(membership_F1_pooled = mean(prim$f1), oracle_pooled = mean(prim$oracle))
  for (b in unique(s$base)) {
    sb <- s[s$base == b, ]; gb <- prim[prim$base == b, ]
    nul <- sb$any_false[sb$scenario == "null"]
    out[[paste0(b, "_F1")]] <- mean(gb$f1)
    out[[paste0(b, "_null_fwer")]] <- mean(nul); out[[paste0(b, "_null_cp_upper")]] <- cp_up(nul)
    out[[paste0(b, "_genes_only_false")]] <- mean(sb$any_false[sb$scenario == "genes_only"])
    out[[paste0(b, "_ambient_decoy")]] <- mean(sb$any_decoy[sb$scenario == "ambient"])
    out[[paste0(b, "_site_any")]] <- mean(sb$any_site, na.rm = TRUE)
    out[[paste0(b, "_within_ctx_false")]] <- mean(sb$any_false[sb$scenario == "context" & sb$mode == "within"])
    out[[paste0(b, "_spurious_pos")]] <- mean(sb$any_spurious[sb$scenario %in% POS_SCEN & sb$mode == "marginal"])
    for (sc in POS_SCEN) out[[paste0(b, "_F1_", sc)]] <- mean(gb$f1[gb$scenario == sc])
  }
  out$median_runtime <- stats::median(s$runtime); out$n_errors <- sum(!is.na(s$error))
  as.data.frame(lapply(out, function(v) round(v, 4)))
}
