## scoring.R -- node-to-truth assignment and cluster scoring for the coordination benchmark.
## Ported from reports/biomni_transfer/part1_callers_scoring.md section 3 (node_assign / eval_calls),
## adapted to sim_coord.R's truth: signals are the planted recurring programs (R100/R50/R20) plus the
## SITE nuisance; private programs are not signals, so a node matching only a private program is
## unassigned and its cluster is spurious.
##
## Node-to-truth rule (PREREG blocker 3): each node (ct, program) is assigned to the signal with the
## highest |cor| between its scores and the signal's true z on the shared donors (member cell types
## only), or to SITE by the multiple R^2 of lm(score ~ strata), whenever the best |cor| >= NODE_THR
## (0.4). A cluster's label is the majority signal over its member nodes (>= 50%), else "spurious";
## a SITE-majority cluster is a "decoy" call.

NODE_THR <- 0.4

## data.frame(cell_type, program_id, signal, cor, type) for every node of the fit
node_assign_coord <- function(fit, truth, strata) {
  out <- list()
  for (ct in fit$cell_types) {
    S <- fit$scores[[ct]]
    if (is.null(S) || !ncol(S)) next
    d <- rownames(S)
    cand <- names(truth$programs)[vapply(truth$programs, function(p) ct %in% p$members, TRUE)]
    C <- matrix(0, ncol(S), length(cand) + 1L, dimnames = list(colnames(S), c(cand, "SITE")))
    for (s in cand) {
      z <- truth$programs[[s]]$z[[ct]][d]
      C[, s] <- abs(suppressWarnings(stats::cor(S, z, use = "pairwise")))
    }
    f <- factor(strata[d])
    C[, "SITE"] <- if (nlevels(f) > 1L) {
      sqrt(pmax(0, apply(S, 2, function(v) summary(stats::lm(v ~ f))$r.squared)))
    } else 0
    C[is.na(C)] <- 0
    best <- colnames(C)[max.col(C, ties.method = "first")]
    bv <- apply(C, 1, max)
    out[[ct]] <- data.frame(cell_type = ct, program_id = colnames(S),
                            signal = ifelse(bv >= NODE_THR, best, NA_character_),
                            cor = bv, stringsAsFactors = FALSE)
  }
  na <- do.call(rbind, out)
  if (is.null(na)) return(na)
  na$type <- ifelse(is.na(na$signal), NA_character_, ifelse(na$signal == "SITE", "site", "target"))
  na
}

## Per planted program: best-cluster F1, plus dataset-level false-call flags.
## Programs with a $type field other than "target" (the artifact decoys "artifact_vdj",
## "artifact_deletion") are excluded from F1/found; clusters whose majority signal is such a
## program are labelled "artifact" (reported like decoys, never scored as targets).
eval_calls_coord <- function(cl, fit, truth, strata) {
  na <- node_assign_coord(fit, truth, strata)
  ptype <- vapply(truth$programs, function(p) p$type %||% "target", "")
  tgt <- names(ptype)[ptype == "target"]
  key <- paste(na$cell_type, na$program_id)
  cmps <- if (nrow(cl)) split(cl, cl$cluster) else list()
  lab <- vapply(cmps, function(c) {
    s <- na$signal[match(paste(c$cell_type, c$program_id), key)]
    tb <- table(s, useNA = "no")
    if (!length(tb) || max(tb) / nrow(c) <= 0.5) "spurious" else names(tb)[which.max(tb)]
  }, "")
  lab_type <- ifelse(lab == "spurious", "spurious",
                     ifelse(lab == "SITE", "decoy",
                            ifelse(ptype[lab] == "target", "target", "artifact")))
  f1 <- found <- stats::setNames(vector("list", length(tgt)), tgt)
  for (s in tgt) {
    mem <- truth$programs[[s]]$members
    per <- vapply(cmps, function(c) {
      s_n <- na$signal[match(paste(c$cell_type, c$program_id), key)]
      ok <- !is.na(s_n) & s_n == s
      rec <- length(unique(c$cell_type[ok])) / length(mem)
      pre <- mean(ok)
      if (rec + pre == 0) 0 else 2 * rec * pre / (rec + pre)
    }, 1)
    f1[[s]] <- if (length(per)) max(per) else 0
    per_found <- vapply(cmps, function(c) {
      s_n <- na$signal[match(paste(c$cell_type, c$program_id), key)]
      ok <- !is.na(s_n) & s_n == s
      length(unique(c$cell_type[ok])) / length(mem) >= 0.5 && mean(ok) >= 0.5 && sum(ok) >= 2L
    }, TRUE)
    found[[s]] <- any(per_found)
  }
  f1 <- unlist(f1); found <- unlist(found)
  brd <- function(sig) {
    if (grepl("R20", sig)) "R20" else if (grepl("R50", sig)) "R50" else "R100"
  }
  list(f1 = f1, found = found, labels = lab, label_type = lab_type,
       summary = data.frame(
         n_cmp = length(cmps),
         n_decoy = sum(lab_type == "decoy"), n_spurious = sum(lab_type == "spurious"),
         n_artifact = sum(lab_type == "artifact"),
         any_decoy = any(lab_type == "decoy"), any_spurious = any(lab_type == "spurious"),
         any_artifact = any(lab_type == "artifact"),
         mean_f1 = if (length(f1)) mean(f1) else NA_real_,
         stringsAsFactors = FALSE))
}
