## summarize.R [results_root] [out_csv]
## Applies the decision rules in PREREG.md to results/<base>/<scenario>/rep*_{ct,pair}.csv.
## Paired over reps (same dataset in every arm); CIs by paired bootstrap over reps.
a <- commandArgs(TRUE)
here <- local({ f <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]); if (is.na(f)) getwd() else dirname(normalizePath(f)) })
res_root <- if (length(a) >= 1) a[1] else file.path(here, "results")
out_csv <- if (length(a) >= 2) a[2] else file.path(res_root, "decisions.csv")
B <- 4000L; set.seed(20260101)

K_ARMS <- c("k_n", "k_20", "k_stab"); R_ARMS <- c("genes_after", "vt_12", "w_S")
MARGIN <- c(P = 0.10, R = 0.02); GUARD <- c(F = 0.03, R = 0.02, time = 3)

rd <- function(base, kind) {
  f <- list.files(file.path(res_root, base), pattern = paste0("_", kind, "\\.csv$"), recursive = TRUE, full.names = TRUE)
  if (!length(f)) return(NULL)
  do.call(rbind, lapply(f, utils::read.csv, stringsAsFactors = FALSE))
}
boot <- function(d) {
  d <- d[is.finite(d)]
  if (length(d) < 2L) return(c(est = NA_real_, lo = NA_real_, hi = NA_real_, n = length(d)))
  m <- replicate(B, mean(sample(d, replace = TRUE)))
  c(est = mean(d), lo = unname(quantile(m, .025)), hi = unname(quantile(m, .975)), n = length(d))
}
TRUE_CTS <- list(S1 = c("A", "B"), S2 = c("A", "C"), S3 = c("C", "D"), S7 = c("A", "B"))

per_rep <- function(base) {
  ct <- rd(base, "ct"); pr <- rd(base, "pair")
  if (is.null(ct) || is.null(pr)) return(NULL)
  chance <- aggregate(chance ~ arm + ct, ct, mean); names(chance)[3] <- "chance_mean"
  ct <- merge(ct, chance, by = c("arm", "ct"))
  ct$Rraw <- ct$best_cor - ct$chance_mean
  R <- do.call(rbind, lapply(split(ct, list(ct$scenario, ct$rep, ct$arm), drop = TRUE), function(g) {
    s <- g$scenario[1]
    if (!s %in% names(TRUE_CTS)) return(NULL)
    gg <- g[g$ct %in% TRUE_CTS[[s]], ]
    data.frame(scenario = s, rep = g$rep[1], arm = g$arm[1], R = mean(gg$Rraw), secs = g$fit_secs[1])
  }))
  P <- aggregate(shared ~ scenario + rep + arm, pr[pr$is_true, ], function(v) as.numeric(any(v)))
  names(P)[4] <- "P"
  Fm <- aggregate(shared ~ scenario + rep + arm, pr[pr$scenario %in% c("S0", "S7"), ], function(v) as.numeric(any(v)))
  names(Fm)[4] <- "F"
  Fo <- aggregate(shared ~ scenario + rep + arm, pr[!pr$is_true & pr$scenario %in% c("S1", "S2", "S3"), ], function(v) as.numeric(any(v)))
  names(Fo)[4] <- "F_other"
  list(R = R, P = P, F = Fm, F_other = Fo, pair = pr)
}
paired <- function(tab, col, arm, scen) {
  a <- tab[tab$arm == arm & tab$scenario %in% scen, c("scenario", "rep", col)]
  b <- tab[tab$arm == "base" & tab$scenario %in% scen, c("scenario", "rep", col)]
  m <- merge(a, b, by = c("scenario", "rep"), suffixes = c("_a", "_b"))
  m[[paste0(col, "_a")]] - m[[paste0(col, "_b")]]
}
mean_of <- function(tab, col, arm, scen) mean(tab[[col]][tab$arm == arm & tab$scenario %in% scen])

rows <- list(); calib <- list()
for (base in c("semi", "param")) {
  pr <- per_rep(base); if (is.null(pr)) next
  for (s in c("S0", "S7")) {
    v <- pr$F[pr$F$arm == "base" & pr$F$scenario == s, "F"]
    if (length(v)) {
      ci <- binom.test(sum(v), length(v))$conf.int
      calib[[length(calib) + 1L]] <- data.frame(base = base, scenario = s, n = length(v), F_base = mean(v), cp_lo = ci[1], cp_hi = ci[2])
    }
  }
  tmed <- tapply(pr$R$secs, pr$R$arm, median)
  for (arm in c(K_ARMS, R_ARMS, "strata")) {
    g_F <- sapply(c("S0", "S7"), function(s) mean_of(pr$F, "F", arm, s) - mean_of(pr$F, "F", "base", s))
    g_R <- if (arm == "strata") 0 else boot(paired(pr$R, "R", arm, c("S1", "S2")))["est"]
    g_t <- unname(tmed[arm] / tmed["base"])
    guards <- is.finite(g_t) && all(g_F <= GUARD["F"], na.rm = TRUE) && g_R >= -GUARD["R"] && g_t <= GUARD["time"]
    if (arm %in% K_ARMS) {
      e <- boot(paired(pr$P, "P", arm, "S3")); thr <- MARGIN["P"]; endpoint <- "dP_S3"
      pass <- is.finite(e["est"]) && e["est"] >= thr && e["lo"] > 0
    } else if (arm %in% R_ARMS) {
      e <- boot(paired(pr$R, "R", arm, c("S1", "S2"))); thr <- MARGIN["R"]; endpoint <- "dR_S1S2"
      pass <- is.finite(e["est"]) && e["est"] >= thr && e["lo"] > 0
    } else {
      e <- boot(paired(pr$F, "F", arm, "S0")); thr <- NA_real_; endpoint <- "dF_S0"
      infl <- mean_of(pr$F, "F", "base", "S0") > 0.05
      cost <- -boot(unlist(lapply(c("S1", "S2", "S3"), function(s) paired(pr$P, "P", arm, s))))["est"]
      pass <- infl && is.finite(e["est"]) && e["hi"] < 0 && is.finite(cost) && cost <= 0.05
    }
    rows[[length(rows) + 1L]] <- data.frame(base = base, arm = arm, endpoint = endpoint, est = e["est"], lo = e["lo"], hi = e["hi"],
      n_reps = e["n"], threshold = thr, effect_pass = pass, guards_pass = guards, dF_S0 = g_F[1], dF_S7 = g_F[2],
      dR = unname(g_R), time_ratio = g_t, row.names = NULL)
  }
}
if (!length(rows)) stop("no results found under ", res_root)
d <- do.call(rbind, rows)
w <- reshape(d[, c("base", "arm", "effect_pass", "guards_pass", "est")], idvar = "arm", timevar = "base", direction = "wide")
sem <- function(col) if (paste0(col, ".semi") %in% names(w)) w[[paste0(col, ".semi")]] else NA
pgetp <- function(col) if (paste0(col, ".param") %in% names(w)) w[[paste0(col, ".param")]] else NA
sign_ok <- ifelse(w$arm == "strata", is.na(pgetp("est")) | pgetp("est") <= 0, pgetp("est") > 0)
w$adopt <- sem("effect_pass") & sem("guards_pass") & pgetp("guards_pass") & sign_ok
w$adopt[is.na(w$adopt)] <- FALSE
w$group <- ifelse(w$arm %in% K_ARMS, "K", ifelse(w$arm %in% c("vt_12", "w_S"), "weights", ifelse(w$arm == "genes_after", "genes", "null")))
order_pref <- c("k_n", "k_20", "k_stab", "vt_12", "w_S", "genes_after", "strata")
w$final <- FALSE
for (gp in unique(w$group)) {
  cand <- w[w$group == gp & w$adopt, "arm"]
  if (length(cand)) w$final[w$arm == cand[order(match(cand, order_pref))][1]] <- TRUE
}
d <- merge(d, w[, c("arm", "adopt", "final")], by = "arm")
dir.create(dirname(out_csv), showWarnings = FALSE, recursive = TRUE)
write.csv(d[order(d$arm, d$base), ], out_csv, row.names = FALSE)
if (length(calib)) write.csv(do.call(rbind, calib), sub("\\.csv$", "_base_calibration.csv", out_csv), row.names = FALSE)
print(w[, c("arm", "adopt", "final")]); if (length(calib)) print(do.call(rbind, calib))
