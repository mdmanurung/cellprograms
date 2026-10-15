suppressPackageStartupMessages(library(data.table))
POS <- c("cmp_distinct", "cmp_mixed", "context", "ambient", "lr", "weak")
brd <- function(sig) fifelse(grepl("R20", sig), "R20", fifelse(grepl("R50", sig), "R50", fifelse(sig == "CTX", "CTX", "R100")))
ld <- function(f) { e <- readRDS(f); g <- as.data.table(e$signals); g[, sr := seed %% 100000]; g[scenario %in% POS & mode == "marginal"] }
G <- rbindlist(list(`cap10` = ld("/mnt/shared-workspace/cmp/evals/heldout2_v2/main/n_agg_noanc.rds"),
                    `cap30 (winner)` = ld("/mnt/shared-workspace/cmp/evals/heldout2_v2/K30/n_agg_noanc.rds"),
                    `node_cor` = ld("/mnt/shared-workspace/cmp/evals/heldout2_v2/main/b_node_cor.rds"),
                    `pa_cc` = ld("/mnt/shared-workspace/cmp/evals/heldout2_v2/main/b_pa_cc.rds")), idcol = "method")
G[, br := brd(signal)]
out <- G[, .(F1 = round(mean(f1), 3), miss = round(mean(f1 == 0), 3), oracle = round(mean(oracle), 3), n = .N),
         by = .(method, base, br)][order(method, base, br)]
print(out, nrows = 100)
