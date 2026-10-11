## prereg_check.R -- evaluate the frozen pre-registered success criteria on the test split.
suppressPackageStartupMessages(library(data.table))
D <- "/mnt/shared-workspace/cmp/evals/test_v11/main"; O <- "/mnt/results/tusoskill/cmp_v1"
cands <- c("n_agg_noanc", "b_node_cor", "b_pa_cc", "e_rf_maxT_gr_amb3", "n_agg_t2", "n_agg_minp_lt35", "b_gene_match", "s_shuffled")
sc <- rbindlist(lapply(cands, function(cn) { f <- file.path(D, paste0(cn, "_score.csv")); if (file.exists(f)) fread(f) }), fill = TRUE)
fwrite(sc, "/tmp/test_sc.csv"); system2("cp", c("/tmp/test_sc.csv", file.path(O, "test_scorecard_v1.1.csv")))
per_ds <- function(cn) { g <- as.data.table(readRDS(file.path(D, paste0(cn, ".rds")))$signals)[mode == "marginal"]
  g[, .(f1 = mean(f1)), by = .(base, scenario, seed)] }
paired <- function(a, b) { m <- merge(per_ds(a), per_ds(b), by = c("base", "scenario", "seed")); d <- m$f1.x - m$f1.y
  set.seed(1); bs <- replicate(5000, mean(sample(d, replace = TRUE)))
  data.table(finalist = a, comparator = b, diff = mean(d), lo = quantile(bs, .025), hi = quantile(bs, .975), n = length(d),
             diff_param = mean(d[m$base == "param"]), diff_semi = mean(d[m$base == "semi"])) }
pr <- rbindlist(lapply(setdiff(cands, "n_agg_noanc"), function(b) paired("n_agg_noanc", b)))
fwrite(pr, "/tmp/test_paired.csv"); system2("cp", c("/tmp/test_paired.csv", file.path(O, "test_paired_vs_finalist.csv")))
f <- sc[candidate == "n_agg_noanc"]; num <- function(x) as.numeric(x)
crit <- data.table(
  criterion = c("1a F1 > b_node_cor (CI>0)", "1b F1 > b_pa_cc (CI>0)", "2 null FWER <= 0.10 (param)", "2 null FWER <= 0.10 (semi)",
                "3 ambient decoy <= 0.10 (param)", "3 ambient decoy <= 0.10 (semi)", "4 within-ctx false <= 0.10 (param)",
                "4 within-ctx false <= 0.10 (semi)", "5 genes_only false <= 0.10 (param)", "5 genes_only false <= 0.10 (semi)",
                "6 s_shuffled F1 < 0.02"),
  value = c(sprintf("%.4f [%.4f, %.4f]", pr[comparator == "b_node_cor"]$diff, pr[comparator == "b_node_cor"]$lo, pr[comparator == "b_node_cor"]$hi),
            sprintf("%.4f [%.4f, %.4f]", pr[comparator == "b_pa_cc"]$diff, pr[comparator == "b_pa_cc"]$lo, pr[comparator == "b_pa_cc"]$hi),
            f$param_null_fwer, f$semi_null_fwer, f$param_ambient_decoy, f$semi_ambient_decoy, f$param_within_ctx_false,
            f$semi_within_ctx_false, f$param_genes_only_false, f$semi_genes_only_false, sc[candidate == "s_shuffled"]$membership_F1_pooled),
  pass = c(pr[comparator == "b_node_cor"]$lo > 0, pr[comparator == "b_pa_cc"]$lo > 0, num(f$param_null_fwer) <= 0.10, num(f$semi_null_fwer) <= 0.10,
           num(f$param_ambient_decoy) <= 0.10, num(f$semi_ambient_decoy) <= 0.10, num(f$param_within_ctx_false) <= 0.10,
           num(f$semi_within_ctx_false) <= 0.10, num(f$param_genes_only_false) <= 0.10, num(f$semi_genes_only_false) <= 0.10,
           num(sc[candidate == "s_shuffled"]$membership_F1_pooled) < 0.02))
fwrite(crit, "/tmp/test_crit.csv"); system2("cp", c("/tmp/test_crit.csv", file.path(O, "test_prereg_criteria.csv")))
print(sc[, .(candidate, F1 = membership_F1_pooled, pF1 = param_F1, sF1 = semi_F1, pFWER = param_null_fwer, sFWER = semi_null_fwer,
             pAmb = param_ambient_decoy, sAmb = semi_ambient_decoy, pWin = param_within_ctx_false, sWin = semi_within_ctx_false,
             pGO = param_genes_only_false, sGO = semi_genes_only_false)])
print(pr); print(crit)
