## scoring.R checks (v2 membership-F1 scorer):
##   R_LIBS_USER=/workspace/.Rlib Rscript benchmarks/coordination/test_scoring.R
here <- local({ f <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]); dirname(normalizePath(f)) })
root <- normalizePath(file.path(here, "..", ".."))
source(file.path(here, "scoring.R"))
if (!requireNamespace("clue", quietly = TRUE)) stop("clue is required for the Hungarian assignment")
chk <- function(ok, msg) if (!isTRUE(ok)) stop("FAILED: ", msg) else cat("ok:", msg, "\n")

## hand-built membership cases: perfect R50, partial R20, spurious extra call
score_selfcheck()

## decoy + trap: a spurious 2-ct call whose donor activity tracks the nuisance, and whose
## members are exactly a priv_corr decoy pair, is both a decoy and a trap
set.seed(7)
bz <- setNames(rnorm(100), sprintf("d%03d", 1:100))
tr <- list(programs = list(list(id = "R50", members = c("ct01", "ct03", "ct05", "ct08", "ct10"))),
           decoy_pairs = list(list(id = "dp1", members = c("ct02", "ct09"))),
           nuisance = bz)
cl <- data.frame(cluster = c(1, 1, 1, 1, 1, 2, 2),
                 cell_type = c("ct01", "ct03", "ct05", "ct08", "ct10", "ct02", "ct09"),
                 program_id = NA_character_, stat = 1, fwer_p = NA, stringsAsFactors = FALSE)
act <- cbind(cluster_1 = setNames(rnorm(100), names(bz)),
             cluster_2 = setNames(bz + rnorm(100, sd = 0.05), names(bz)))
s <- score_method(list(method = "chk", clusters = cl, activity = act), tr)
chk(isTRUE(s$program_level$found[s$program_level$program_id == "R50"]), "perfect R50 call found")
chk(s$dataset_level$n_calls == 2L && s$dataset_level$n_spurious == 1L, "one accepted call, one spurious")
chk(s$dataset_level$n_decoy == 1L, "spurious call with nuisance-correlated activity is a decoy")
chk(s$dataset_level$n_trap == 1L, "spurious call containing both priv_corr pair members is a trap")

## same spurious call without nuisance-correlated activity: trap fires, decoy does not
act0 <- cbind(cluster_1 = setNames(rnorm(100), names(bz)),
              cluster_2 = setNames(rnorm(100), names(bz)))
s0 <- score_method(list(method = "chk", clusters = cl, activity = act0), tr)
chk(s0$dataset_level$n_decoy == 0L && s0$dataset_level$n_trap == 1L, "uncorrelated activity: no decoy, trap still fires")

## empty call set: no calls, the planted program is scored not-found
cl0 <- data.frame(cluster = integer(), cell_type = character(), program_id = character(),
                  stat = numeric(), fwer_p = numeric(), stringsAsFactors = FALSE)
e <- score_method(list(method = "chk", clusters = cl0, activity = NULL), tr)
chk(e$dataset_level$n_calls == 0L && nrow(e$program_level) == 1L && !e$program_level$found[[1]],
    "empty call set: program scored not found")

cat("all scoring checks passed\n")
