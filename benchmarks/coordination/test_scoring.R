## scoring.R checks: R_LIBS_USER=/nonexistent Rscript benchmarks/coordination/test_scoring.R
here <- local({ f <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]); dirname(normalizePath(f)) })
root <- normalizePath(file.path(here, "..", ".."))
library(stats); library(utils)
for (f in list.files(file.path(root, "R"), full.names = TRUE)) source(f)
source(file.path(here, "methods.R"))
source(file.path(here, "scoring.R"))
set.seed(3)
cts <- sprintf("ct%02d", 1:10); ids <- sprintf("d%03d", 1:100)
strata <- setNames(sample(c("a", "b", "c"), 100, TRUE), ids)
bz <- setNames(c(a = 0, b = 1.5, c = 3)[strata] + rnorm(100, sd = 0.1), ids)  # nuisance = site-driven
mem <- list(R100 = cts, R50 = cts[1:5], R20 = cts[c(9, 10)])
zl <- lapply(mem, function(m) setNames(rnorm(100), ids))
sc <- list()
for (ct in cts) {
  d <- ids
  Z <- matrix(rnorm(100 * 4), 100, 4, dimnames = list(d, paste0(ct, "_", 1:4)))
  for (s in seq_along(names(mem))) if (ct %in% mem[[s]]) Z[, s] <- 2 * zl[[s]][d] + rnorm(100, sd = 0.5)
  if (ct == "ct01") Z[, 4] <- 3 * bz[d]              # pure nuisance node -> SITE/decoy
  sc[[ct]] <- Z
}
fit <- list(cell_types = cts, scores = sc)
truth <- list(programs = setNames(lapply(names(mem), function(s) {
  list(members = mem[[s]], z = setNames(lapply(mem[[s]], function(ct) zl[[s]]), mem[[s]]))
}), names(mem)))

na <- node_assign_coord(fit, truth, strata)
chk <- function(ok, msg) if (!isTRUE(ok)) stop("FAILED: ", msg) else cat("ok:", msg, "\n")
chk(all(na$signal[na$cell_type == "ct03" & na$program_id == "ct03_1"] == "R100"), "planted node assigned to its signal")
chk(all(na$signal[na$cell_type == "ct03" & na$program_id == "ct03_4"] == "R50" |
        is.na(na$signal[na$cell_type == "ct03" & na$program_id == "ct03_4"])), "non-signal node unassigned or weak")
ct01_site <- na[na$cell_type == "ct01" & na$program_id == "ct01_4", ]
chk(identical(ct01_site$signal, "SITE") && identical(ct01_site$type, "site"), "nuisance node assigned to SITE")
## non-member nodes must not match a signal they are not part of
chk(all(is.na(na$signal[na$cell_type == "ct06" & na$program_id %in% c("ct06_2", "ct06_3")])), "non-member nodes unassigned")

## eval_calls: a cluster spanning the R50 members scores F1 = 1 and found = TRUE
cl <- data.frame(cluster = 1, cell_type = mem$R50,
                 program_id = paste0(mem$R50, "_2"), stat = 5, fwer_p = 0.01, stringsAsFactors = FALSE)
ev <- eval_calls_coord(cl, fit, truth, strata)
chk(abs(ev$f1[["R50"]] - 1) < 1e-12, "perfect R50 cluster: F1 = 1")
chk(isTRUE(ev$found[["R50"]]), "perfect R50 cluster: found")
chk(ev$summary$n_decoy == 0, "no decoy calls for a true cluster")
## the nuisance-only cluster is a decoy
cld <- data.frame(cluster = 1, cell_type = "ct01", program_id = "ct01_4", stat = 5, fwer_p = 0.01, stringsAsFactors = FALSE)
evd <- eval_calls_coord(cld, fit, truth, strata)
chk(evd$summary$n_decoy == 1 && isTRUE(evd$summary$any_decoy), "SITE-majority cluster labelled decoy")
chk(all(ev$f1 == 0) || !isTRUE(ev$found[["R20"]]), "no R20 cluster formed")
## empty call set: F1 = 0, nothing found
eve <- eval_calls_coord(.empty_clusters(), fit, truth, strata)
chk(all(eve$f1 == 0) && eve$summary$n_cmp == 0L, "empty clusters: F1 = 0")
cat("all scoring checks passed\n")
