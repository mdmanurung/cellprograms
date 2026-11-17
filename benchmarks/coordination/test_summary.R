## test_summary.R -- add_family() and an end-to-end run of summarize_main.R on synthetic CSVs.
## run: Rscript benchmarks/coordination/test_summary.R
here <- local({ f <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]); dirname(normalizePath(f)) })
source(file.path(here, "family.R"))
chk <- function(ok, msg) if (!isTRUE(ok)) stop("FAILED: ", msg) else cat("ok:", msg, "\n")

## add_family on a hand-made table
pd <- data.frame(method = c("a", "b"), n_pairs = c(64L, 64L), mean_diff = c(0.1, -0.05),
                 ci_lo = c(0.02, -0.15), ci_hi = c(0.18, 0.05))
fam <- add_family(pd, "demo")
chk(identical(fam$family_size, c(2L, 2L)) && all(fam$family == "demo") && all(fam$adjustment == "none"),
    "family columns name the family, its size and 'none' adjustment")
se <- (0.18 - 0.1) / 1.96
chk(abs(fam$p_unadjusted[1] - 2 * pt(-0.1 / se, 63)) < 1e-12, "p_unadjusted is the two-sided t p-value from mean and CI")
chk(fam$p_unadjusted[1] < 0.05 && fam$p_unadjusted[2] > 0.05, "p-values order the two contrasts")
chk(nrow(add_family(pd[0, ], "demo")) == 0, "empty table passes through")

## end to end: copy the script next to a synthetic out/ tree and run it
tmp <- file.path(tempdir(), "sm"); dir.create(file.path(tmp, "out", "t", "prog"), recursive = TRUE)
dir.create(file.path(tmp, "out", "t", "ds"))
file.copy(file.path(here, c("summarize_main.R", "family.R")), tmp)
set.seed(1)
prog <- do.call(rbind, lapply(c("pa_cc", "node_cor"), function(m) data.frame(
  mode = "act_only", breadth = rep(c("R50", "R20"), each = 20), amp = 0.5, act_cor = 1, rep = rep(1:20, 2),
  program_id = rep(c("R50", "R20"), each = 20), method = m, cell_id = "act_only|x",
  found = TRUE, recall = pmin(1, runif(40, 0.3, 0.9) + (m == "node_cor") * 0.05), precision = 1, breadth_error = 0)))
ds <- do.call(rbind, lapply(c("pa_cc", "node_cor"), function(m) data.frame(
  mode = "act_only", breadth = "R50", amp = 0.5, act_cor = 1, method = m, n_spurious = 0L,
  spurious_per_call = 0, n_decoy = 0L, n_trap = 0L, runtime_secs = 1)))
write.csv(prog, file.path(tmp, "out", "t", "prog", "a.csv"), row.names = FALSE)
write.csv(ds, file.path(tmp, "out", "t", "ds", "a.csv"), row.names = FALSE)
out <- system2(file.path(R.home("bin"), "Rscript"), c(shQuote(file.path(tmp, "summarize_main.R"))),
               env = "OUT_SUBDIR=t", stdout = TRUE, stderr = TRUE)
res <- read.csv(file.path(tmp, "out", "t", "summary_paired_primary.csv"))
chk(all(c("p_unadjusted", "family", "family_size", "adjustment") %in% names(res)), "summarize_main.R writes the family columns")
chk(res$family_size == 1L && res$method == "node_cor", "family of one contrast: node_cor vs pa_cc")
chk(any(grepl("1 arm-versus-pa_cc contrasts, p-values unadjusted", out, fixed = TRUE)), "the printed table names the family size")
cat("all summary checks passed\n")
