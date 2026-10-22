## test_main_run.R — smoke test for the replicate harness (fast methods only).
## run: R_LIBS_USER=/workspace/.Rlib Rscript benchmarks/coordination/test_main_run.R
here <- local({ f <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]); dirname(normalizePath(f)) })
root <- normalizePath(file.path(here, "..", ".."))
library(stats)
for (f in c("data-model.R", "utils.R", "ebmf.R", "principal-angles.R", "block-joint.R", "eb-cp.R")) source(file.path(root, "R", f))
for (f in c("sim_coord.R", "methods.R", "eb_cp.R", "scoring.R", "replicate.R")) source(file.path(here, f))
chk <- function(ok, msg) if (!isTRUE(ok)) stop("FAILED: ", msg) else cat("ok:", msg, "\n")

outdir <- file.path(here, "out", "smoke")
dir.create(file.path(outdir, "prog"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(outdir, "ds"), recursive = TRUE, showWarnings = FALSE)
fm <- c("node_cor", "gene_match", "eb_cp")

## null replicate: ds rows for all methods, header-only prog CSV
run_one_replicate("null", "R50", 0.5, 1L, 4242L, 1, 2000L, outdir, methods = fm)
ds <- read.csv(file.path(outdir, "ds", "null_R50_a0.5_c1_r001.csv"))
chk(nrow(ds) == length(fm), "null: one ds row per method")
chk(all(ds$n_calls == 0), "null: zero calls from all fast methods")
pr <- read.csv(file.path(outdir, "prog", "null_R50_a0.5_c1_r001.csv"))
chk(nrow(pr) == 0L && "program_id" %in% names(pr), "null: header-only prog CSV")

## signal replicate: eb_cp should find the planted R100 (act+genes regime is
## gate-test-proven); prog CSV populated
run_one_replicate("act+genes", "R100", 0.75, 1L, 4343L, 1, 2000L, outdir, methods = fm)
ds2 <- read.csv(file.path(outdir, "ds", "act+genes_R100_a0.75_c1_r001.csv"))
pr2 <- read.csv(file.path(outdir, "prog", "act+genes_R100_a0.75_c1_r001.csv"))
chk(nrow(ds2) == length(fm), "signal: one ds row per method")
chk(nrow(pr2) == length(fm), "signal: one prog row per method (one planted program)")
chk(any(pr2$method == "eb_cp" & pr2$found), "signal: eb_cp finds R100 at amp 0.75")

## resumability: second call skips
st <- run_one_replicate("null", "R50", 0.5, 1L, 4242L, 1, 2000L, outdir, methods = fm)
chk(identical(st, "skip") || identical(as.character(st), "skip"), "resumability: existing replicate skipped")

## error isolation: force a method failure via a bogus method name in a copy
## (run_method stops on unknown method -> tryCatch records an error row)
run_one_replicate("null", "R50", 0.5, 2L, 4244L, 1, 2000L, outdir,
                  methods = c("node_cor", "bogus_method"))
ds3 <- read.csv(file.path(outdir, "ds", "null_R50_a0.5_c1_r002.csv"))
chk(nrow(ds3) == 2L, "error isolation: both methods produce ds rows")
chk(sum(!is.na(ds3$error)) == 1L, "error isolation: exactly one error row")
chk(all(c("n_calls", "n_spurious", "stage1_secs") %in% names(ds3)),
    "error isolation: error rows carry the full column set")

cat("main-run harness smoke checks passed\n")
