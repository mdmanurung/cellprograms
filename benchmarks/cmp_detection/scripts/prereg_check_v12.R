## prereg_check_v12.R -- evaluate the frozen cmp-bench-v1.2 criteria on held-out seeds 1051-1100.
## Winner: n_agg_noanc on K30 fits. Comparators: n_agg_noanc / b_node_cor / b_pa_cc on main (cap-10) fits.
## Sanity: s_shuffled on K30 fits. Paired bootstrap: 2000 resamples (as pre-registered), seed 1.
suppressPackageStartupMessages(library(data.table))
E <- "/mnt/shared-workspace/cmp/evals/heldout2_v2"; O <- "/mnt/results/tusoskill/cmp_v2"
wr <- function(x, f) { tmp <- tempfile(fileext = ".csv"); fwrite(x, tmp); system2("cp", c(tmp, file.path(O, f))) }
runs <- list(winner_K30 = c("K30", "n_agg_noanc"), finalist_K10 = c("main", "n_agg_noanc"),
             b_node_cor = c("main", "b_node_cor"), b_pa_cc = c("main", "b_pa_cc"), s_shuffled_K30 = c("K30", "s_shuffled"))
sc <- rbindlist(lapply(names(runs), function(k) cbind(label = k, fread(file.path(E, runs[[k]][1], paste0(runs[[k]][2], "_score.csv"))))), fill = TRUE)
wr(sc, "heldout_v1.2_scorecard.csv")
per_ds <- function(k) { g <- as.data.table(readRDS(file.path(E, runs[[k]][1], paste0(runs[[k]][2], ".rds")))$signals)[mode == "marginal"]
  g[, .(f1 = mean(f1)), by = .(base, scenario, seed)] }
W <- per_ds("winner_K30")
paired <- function(b) { m <- merge(W, per_ds(b), by = c("base", "scenario", "seed")); d <- m$f1.x - m$f1.y
  set.seed(1); bs <- replicate(2000, mean(sample(d, replace = TRUE)))
  data.table(winner = "K30+n_agg_noanc", comparator = b, diff = mean(d), lo = unname(quantile(bs, .025)), hi = unname(quantile(bs, .975)),
             n = length(d), diff_param = mean(d[m$base == "param"]), diff_semi = mean(d[m$base == "semi"])) }
pr <- rbindlist(lapply(c("finalist_K10", "b_node_cor", "b_pa_cc"), paired))
wr(pr, "heldout_v1.2_paired.csv")
f <- sc[label == "winner_K30"]; num <- function(x) as.numeric(x); ci <- function(b) sprintf("%.4f [%.4f, %.4f]", pr[comparator == b]$diff, pr[comparator == b]$lo, pr[comparator == b]$hi)
crit <- data.table(
  criterion = c("1 F1 > v1 finalist cap10 (CI>0)", "2a F1 > b_node_cor (CI>0)", "2b F1 > b_pa_cc (CI>0)",
                "3 null FWER <= 0.10 (param)", "3 null FWER <= 0.10 (semi)", "4 ambient decoy <= 0.10 (param)", "4 ambient decoy <= 0.10 (semi)",
                "5 within-ctx false <= 0.10 (param)", "5 within-ctx false <= 0.10 (semi)",
                "6 genes_only false <= 0.10 (param)", "6 genes_only false <= 0.10 (semi)", "7 s_shuffled F1 < 0.02"),
  value = c(ci("finalist_K10"), ci("b_node_cor"), ci("b_pa_cc"), f$param_null_fwer, f$semi_null_fwer, f$param_ambient_decoy, f$semi_ambient_decoy,
            f$param_within_ctx_false, f$semi_within_ctx_false, f$param_genes_only_false, f$semi_genes_only_false,
            sc[label == "s_shuffled_K30"]$membership_F1_pooled),
  pass = c(pr[comparator == "finalist_K10"]$lo > 0, pr[comparator == "b_node_cor"]$lo > 0, pr[comparator == "b_pa_cc"]$lo > 0,
           num(f$param_null_fwer) <= 0.10, num(f$semi_null_fwer) <= 0.10, num(f$param_ambient_decoy) <= 0.10, num(f$semi_ambient_decoy) <= 0.10,
           num(f$param_within_ctx_false) <= 0.10, num(f$semi_within_ctx_false) <= 0.10,
           num(f$param_genes_only_false) <= 0.10, num(f$semi_genes_only_false) <= 0.10,
           num(sc[label == "s_shuffled_K30"]$membership_F1_pooled) < 0.02))
wr(crit, "heldout_v1.2_criteria.csv")
cat("n_files per run:\n"); print(sc[, .(label, n_files)])
print(sc[, .(label, F1 = membership_F1_pooled, pF1 = param_F1, sF1 = semi_F1, pFWER = param_null_fwer, sFWER = semi_null_fwer,
             pAmb = param_ambient_decoy, sAmb = semi_ambient_decoy, pWin = param_within_ctx_false, sWin = semi_within_ctx_false,
             pGO = param_genes_only_false, sGO = semi_genes_only_false)])
print(pr); print(crit); cat("ALL PASS:", all(crit$pass), "\n")
