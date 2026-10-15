## quick within-mode re-score on context-scenario fits only (benchmark v1.1 strata = site x context)
args <- commandArgs(TRUE); cands <- strsplit(args[1], ",")[[1]]; cores <- as.integer(args[2])
source("/mnt/shared-workspace/cmp/bench/common.R"); source("/mnt/shared-workspace/cmp/bench/eval_cmp.R"); source("/mnt/shared-workspace/cmp/bench/registry.R")
files <- unlist(lapply(c("param", "semi"), function(b) list.files(file.path("/mnt/shared-workspace/cmp/fits/main", b, "context"), full.names = TRUE)))
for (cn in cands) {
  ev <- eval_method(REG[[cn]], files, cores, cn); s <- as.data.frame(ev$summary); g <- as.data.frame(ev$signals)
  out <- do.call(rbind, lapply(c("param", "semi"), function(b) data.frame(cand = cn, base = b,
    within_false = mean(s$any_false[s$base == b & s$mode == "within"] & s$n_cmp[s$base == b & s$mode == "within"] > 0),
    within_spur = mean(s$any_spurious[s$base == b & s$mode == "within"]),
    within_R20_F1 = mean(g$f1[g$base == b & g$mode == "within" & g$signal == "R20"]),
    marg_F1 = mean(g$f1[g$base == b & g$mode == "marginal"]))))
  print(out); write.csv(out, sprintf("/mnt/shared-workspace/cmp/evals/ctx_test_%s.csv", cn), row.names = FALSE)
}
