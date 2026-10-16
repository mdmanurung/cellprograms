# Aggregate small-batch benchmark results (list of run records) and apply the rule.
suppressPackageStartupMessages(library(cellprograms))

scen <- list(
  list(id = "S1 (32 donors, 8/batch, sigma=0.375)", path = "/workspace/batch_small_s0375/results.rds", gate = TRUE),
  list(id = "S1b (32 donors, 8/batch, sigma=0.75)", path = "/workspace/batch_small_s075/results.rds", gate = FALSE),
  list(id = "S2 (20 donors, 5/batch, sigma=0.375)", path = "/workspace/batch_small_s05/results.rds", gate = FALSE)
)

for (sc in scen) {
  r <- readRDS(sc$path)
  cat("\n===", sc$id, "| runs:", length(r), "===\n")
  arms <- names(r[[1]])
  rec <- list()
  for (a in intersect(c("A0","A1","A2","A3"), arms)) {
    v <- do.call(rbind, lapply(r, function(run) data.frame(
      recovery = run[[a]]$recovery, batch_r2 = run[[a]]$batch_r2, K = run[[a]]$K)))
    rec[[a]] <- data.frame(arm = a,
                           recovery = round(mean(v$recovery), 4), rec_sd = round(sd(v$recovery), 4),
                           batch_r2 = round(mean(v$batch_r2), 5), r2_sd = round(sd(v$batch_r2), 5),
                           K = paste(sort(unique(v$K)), collapse = ","))
  }
  summ <- do.call(rbind, rec)
  print(summ, row.names = FALSE)
  if (sc$gate && all(c("A0", "A1") %in% arms)) {
    v0 <- do.call(rbind, lapply(r, function(run) data.frame(rec = run$A0$recovery, r2 = run$A0$batch_r2, K = run$A0$K)))
    v1 <- do.call(rbind, lapply(r, function(run) data.frame(rec = run$A1$recovery, r2 = run$A1$batch_r2, K = run$A1$K)))
    d_rec <- mean(v1$rec) - mean(v0$rec)
    d_r2 <- mean(v1$r2) - mean(v0$r2)
    k_ok <- all(v1$K <= 2 * v0$K)
    fire <- (d_rec >= 0.02) && (d_r2 <= 0) && k_ok
    cat("S1 gate: recovery delta =", round(d_rec, 4), "| batchR2 delta =", round(d_r2, 5),
        "| K not inflated:", k_ok, "| RULE FIRES:", fire, "\n")
  }
}
