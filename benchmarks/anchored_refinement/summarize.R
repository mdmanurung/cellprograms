## summarize.R -- paired summary + registered decision rule
## (PREREG_mr4-anchored-v1.md): adopt anchored refinement iff, vs denovo_cell,
## (1) F1 diff >= +0.02, (2) paired bootstrap 95% CI excludes 0, (3) no
## reconstruction degradation > 1% (recon_ratio_min <= 1.01).
##   Rscript summarize.R out/
`%||%` <- function(a, b) if (is.null(a)) b else a
arg1 <- commandArgs(TRUE)[1] %||% "out"
d <- do.call(rbind, lapply(list.files(arg1, pattern = "^rep[0-9]+\\.csv$",
                                      full.names = TRUE), read.csv))
wide <- reshape(d[, c("rep", "arm", "f1", "recon_ratio_min")], idvar = "rep",
                timevar = "arm", direction = "wide")
fnames <- grep("^f1\\.", names(wide), value = TRUE)
arms <- sub("^f1\\.", "", fnames)
ref <- "f1.denovo_cell"
boot <- function(x, n = 2000, seed = 1) {
  set.seed(seed); v <- replicate(n, mean(sample(x, replace = TRUE), na.rm = TRUE))
  quantile(v, c(0.025, 0.975), na.rm = TRUE)
}
out <- do.call(rbind, lapply(setdiff(arms, "denovo_cell"), function(a) {
  diff <- wide[[paste0("f1.", a)]] - wide[[ref]]
  ci <- boot(diff)
  rr <- max(wide[[paste0("recon_ratio_min.", a)]], na.rm = TRUE)
  data.frame(arm = a, mean_f1 = mean(wide[[paste0("f1.", a)]], na.rm = TRUE),
             mean_f1_ref = mean(wide[[ref]], na.rm = TRUE),
             d_f1 = mean(diff, na.rm = TRUE), ci_lo = ci[1], ci_hi = ci[2],
             recon_ratio_max = rr,
             adopt = mean(diff, na.rm = TRUE) >= 0.02 && (ci[1] > 0) && rr <= 1.01)
}))
cat("Registered decision (PREREG_mr4-anchored-v1):\n")
print(out, row.names = FALSE)
write.csv(out, file.path(arg1, "decision.csv"), row.names = FALSE)
