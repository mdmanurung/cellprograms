## 06_real_pa_vs_pmd.R -- PA-scores vs multi-view PMD on REAL COMBAT (repo K=10 fit).
##
## Arms: PA-scores (all 45 CT pairs, n_perm=200) and PMD (n_states=5, n_perm=100, seed 43),
## each on raw scores and on scores residualized on sex + Institute (per cell type, lm
## residuals on that cell type's own observed donors; donors with NA covariates keep their
## centred raw score). Source is deliberately NOT residualized (biology of interest).
## Pair-calling rules copied from 05_sim_pa_vs_pmd.R:
##   PA  : shared iff principal_angles()$n_shared_05 >= 1
##   PMD : shared iff some MCP has perm p < 0.05 AND pair_cors[a,b] > 0.5
## Usage: Rscript 06_real_pa_vs_pmd.R   (SMOKE=1 env var: n_perm 10/10, 2 states)

here <- dirname(normalizePath(sub("--file=", "", grep("--file=", commandArgs(), value = TRUE)[1])))
suppressPackageStartupMessages(pkgload::load_all(file.path(here, "fork"), quiet = TRUE))
smoke <- nzchar(Sys.getenv("SMOKE"))
SEED <- 43L; N_PERM_PA <- if (smoke) 10L else 200L; N_PERM_PMD <- if (smoke) 10L else 100L
N_STATES <- if (smoke) 2L else 5L
out <- file.path(here, "results", if (smoke) "real_smoke" else "real"); dir.create(out, FALSE, TRUE)

fit0 <- readRDS(file.path(here, "results/fits/combat_repo_K10.rds"))
meta <- read.csv(file.path(here, "data/combat/donor_meta.csv"), row.names = 1)
age_map <- c("19-30" = 24.5, "31-40" = 35.5, "41-50" = 45.5, "51-60" = 55.5,
             "61-70" = 65.5, "71-80" = 75.5, "81-90" = 85.5, ">=91" = 92.5)
meta$Age_num <- unname(age_map[as.character(meta$Age)])
var <- read.csv(file.path(here, "data/combat/pb_var.csv"))
id2sym <- setNames(var$feature_name, var$feature_id)

## ---- residualization (sex + Institute) ---------------------------------------
residualize <- function(fit) {
  for (ct in names(fit$scores)) {
    Z <- fit$scores[[ct]]
    if (ncol(Z) == 0L) next
    d <- meta[rownames(Z), c("sex", "Institute"), drop = FALSE]
    cc <- stats::complete.cases(d)
    R <- sweep(Z, 2L, colMeans(Z), "-")                     # fallback: centred raw
    d <- droplevels(data.frame(lapply(d[cc, ], factor)))
    terms <- names(Filter(function(f) nlevels(f) > 1L, d))  # drop constant covariates
    if (length(terms) > 0L)
      R[cc, ] <- stats::resid(stats::lm(Z[cc, , drop = FALSE] ~ ., data = d[terms]))
    else R[cc, ] <- sweep(Z[cc, , drop = FALSE], 2L, colMeans(Z[cc, , drop = FALSE]), "-")
    fit$scores[[ct]] <- R
  }
  fit
}

## ---- PA arm ---------------------------------------------------------------------
run_pa <- function(fit) {
  cts <- names(fit$scores)[vapply(fit$scores, ncol, 1L) > 0L]
  pr <- t(utils::combn(cts, 2L))
  do.call(rbind, lapply(seq_len(nrow(pr)), function(i) {
    a <- pr[i, 1]; b <- pr[i, 2]
    pa <- suppressWarnings(principal_angles(fit, a, b, space = "scores",
                                            n_perm = N_PERM_PA, seed = SEED))
    data.frame(pair = paste(a, b, sep = "-"), ct_a = a, ct_b = b,
               k_a = unname(pa$dims["k_a"]), k_b = unname(pa$dims["k_b"]),
               n_donors = unname(pa$dims["effective_dim"]),
               n_shared = pa$n_shared_05,
               min_angle_deg = min(pa$angles) * 180 / pi,
               max_cosine = max(pa$cosines), asymptotic_ref = pa$asymptotic_ref,
               min_p = min(pa$p_values), shared = pa$n_shared_05 >= 1L)
  }))
}

## ---- PMD arm -----------------------------------------------------------------------
## trace_state() has no 'pmd' branch (fork/R/integrate.R), so reuse it through an
## ICA-shaped integ holding the MCP program weights. PMD weights act on STANDARDIZED
## scores, so divide by each program's score SD (own donors) to get gene-space weights.
SIGS <- list(
  Ig        = "^IG[HKL][VCJ]|^JCHAIN$|^MZB1$|^IGLL5$|^TNFRSF17$",
  sex       = "^XIST$|^TSIX$|^RPS4Y1$|^DDX3Y$|^UTY$|^KDM5D$|^EIF1AY$|^ZFY$|^TXLNGY$|^USP9Y$|^NLGN4Y$",
  immed_early = "^FOS$|^FOSB$|^JUN$|^JUNB$|^JUND$|^DUSP1$|^IER2$|^IER3$|^ZFP36$|^NR4A1$|^EGR1$|^KLF6$|^ATF3$|^BTG2$",
  GSTM1     = "^GSTM1$|^GSTM2$|^GSTM4$",
  interferon = "^IFIT[0-9]$|^IFI44L?$|^RSAD2$|^MX[12]$|^ISG15$|^OAS[123L]$|^IFI6$|^XAF1$|^IFI27$|^IFITM3$|^LY6E$|^EPSTI1$")

run_pmd <- function(fit) {
  fit <- suppressWarnings(integrate_programs(fit, method = "pmd", n_states = N_STATES,
                          n_perm = N_PERM_PMD, seed = SEED, verbose = 1L))
  fit
}

pmd_tables <- function(fit, tag) {
  integ <- fit$integration; mcps <- integ$mcps
  sc <- fit$scores
  sdv <- unlist(unname(lapply(sc, function(S) apply(S, 2L, stats::sd))))
  pw <- integ$program_weights
  rows <- list(); assoc <- list(); sigr <- list()
  for (m in mcps) {
    id <- m$state_id; st <- m$state
    pc <- m$pair_cors; iu <- upper.tri(pc)
    contr <- sort(m$view_contribution, decreasing = TRUE)
    top3 <- names(contr)[1:min(3, length(contr))]
    w <- pw[, id]; w <- w / sdv[names(w)]; w[!is.finite(w)] <- 0
    tr30 <- trace_state(fit, integ = list(method = "ica", program_weights = cbind(w) |>
                          `colnames<-`(id)), state = id, top_n = 30L)
    sym <- function(g) id2sym[names(g)]
    fmt <- function(ct, n) { g <- tr30[[ct]]; if (is.null(g)) return(NA_character_)
      g <- g[seq_len(min(n, length(g)))]
      paste0(ifelse(g > 0, "+", "-"), sym(g), collapse = " ") }
    rows[[id]] <- data.frame(arm = tag, state = id, objective = m$objective,
      perm_p = m$pvalue, rho = m$rho, n_active = length(m$active_celltypes),
      mean_pair_r = mean(pc[iu][is.finite(pc[iu])]),
      mean_view_contribution = mean(m$view_contribution[m$active_celltypes]),
      top3_cts = paste(top3, collapse = "+"),
      top8_genes = paste(vapply(top3, function(ct) paste0(ct, ": ", fmt(ct, 8)), ""), collapse = " | "))
    ## associations (donor-level state = mean of per-view MCP scores)
    d <- meta[names(st), ]; ok <- !is.na(st)
    kw <- function(v) { v <- droplevels(factor(v[ok])); v <- v; keep <- !is.na(v)
      if (nlevels(v) < 2L) NA_real_ else stats::kruskal.test(st[ok][keep], v[keep])$p.value }
    sx <- d$sex[ok]
    assoc[[id]] <- data.frame(arm = tag, state = id,
      sex_wilcox_p = stats::wilcox.test(st[ok] ~ factor(sx))$p.value,
      Institute_kruskal_p = kw(d$Institute), Source_kruskal_p = kw(d$Source),
      Outcome_kruskal_p = kw(d$Outcome),
      Age_spearman_rho = unname(stats::cor(st[ok], d$Age_num[ok], method = "spearman", use = "complete.obs")),
      Age_spearman_p = suppressWarnings(stats::cor.test(st[ok], d$Age_num[ok], method = "spearman")$p.value))
    ## signature hits among the top-30 genes of each active CT
    hits <- vapply(names(SIGS), function(s) {
      h <- vapply(tr30, function(g) sum(grepl(SIGS[[s]], sym(g))), 0L)
      sprintf("%d(max/CT)|%dCTs", max(h), sum(h >= 2L)) }, "")
    sigr[[id]] <- data.frame(arm = tag, state = id, t(hits))
  }
  list(mcp = do.call(rbind, rows), assoc = do.call(rbind, assoc), sig = do.call(rbind, sigr))
}

pmd_pairs <- function(fit, tag, pairs) {  # PMD pair call, same rule as 05
  mcps <- fit$integration$mcps
  called <- vapply(strsplit(pairs, "-"), function(p) {
    any(vapply(mcps, function(m) {
      pc <- if (all(p %in% rownames(m$pair_cors))) m$pair_cors[p[1], p[2]] else NA_real_
      isTRUE(m$pvalue < 0.05) && isTRUE(pc > 0.5)
    }, TRUE))
  }, TRUE)
  data.frame(pair = pairs, shared = called, arm = tag)
}

## ---- main --------------------------------------------------------------------------
stopifnot(all(rownames(fit0$scores[[1]]) %in% rownames(meta)))
fits <- list(raw = fit0, resid = residualize(fit0))
res <- list()
for (tag in names(fits)) {
  t0 <- proc.time()[["elapsed"]]
  pa <- run_pa(fits[[tag]])
  pa$shared_count <- sum(pa$shared)
  utils::write.csv(pa, file.path(out, sprintf("pa_pairs_%s.csv", tag)), row.names = FALSE)
  message(sprintf("[%s] PA: %d/%d pairs shared (%.0fs)", tag, sum(pa$shared), nrow(pa),
                  proc.time()[["elapsed"]] - t0))
  t0 <- proc.time()[["elapsed"]]
  f <- run_pmd(fits[[tag]]); tabs <- pmd_tables(f, tag)
  utils::write.csv(tabs$mcp, file.path(out, sprintf("pmd_mcps_%s.csv", tag)), row.names = FALSE)
  utils::write.csv(tabs$assoc, file.path(out, sprintf("pmd_mcp_assoc_%s.csv", tag)), row.names = FALSE)
  utils::write.csv(tabs$sig, file.path(out, sprintf("pmd_mcp_signature_hits_%s.csv", tag)), row.names = FALSE)
  pp <- pmd_pairs(f, tag, pa$pair)
  pp$pa_shared <- pa$shared
  utils::write.csv(pp, file.path(out, sprintf("pmd_pairs_%s.csv", tag)), row.names = FALSE)
  saveRDS(f$integration, file.path(out, sprintf("pmd_integration_%s.rds", tag)))
  message(sprintf("[%s] PMD: %d MCPs, %d/%d pairs shared (%.0fs)", tag, length(f$integration$mcps),
                  sum(pp$shared), nrow(pp), proc.time()[["elapsed"]] - t0))
}
## raw-vs-resid association table side by side
a <- do.call(rbind, lapply(names(fits), function(tg) utils::read.csv(file.path(out, sprintf("pmd_mcp_assoc_%s.csv", tg)))))
utils::write.csv(a, file.path(out, "pmd_mcp_assoc_raw_vs_resid.csv"), row.names = FALSE)
message("done")
