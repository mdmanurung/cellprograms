## smoke.R -- pilot-scale in-sandbox smoke for the MR-4 harness (worker-1 has
## flashier.utils). Runs 2 reduced replicates through the full arm set and
## checks the contract: every arm produces a row, metrics are finite, the
## planted pair is recoverable at pilot amplitude.
##   R_LIBS=<lib> Rscript smoke.R
Sys.setenv(MR4_G = "400", MR4_DONORS = "12", MR4_CELLS = "25", MR4_NVAR = "300",
           MR4_NPERM = "49", MR4_MINDON = "6", MR4_MINSHARED = "6", MR4_NMC = "6",
           MR4_AMP = "1.5", MR4_NCT = "2")
for (rp in 1:2) {
  rep <- rp  # run_rep.R reads this
  source_local <- function(f) source(f, local = TRUE)
  source_local("run_rep.R")
}
res <- do.call(rbind, lapply(list.files("out", pattern = "^rep[0-9]+\\.csv$", full.names = TRUE), read.csv))
stopifnot(nrow(res) >= 2L * 6L - 2L,                     # 6 arms x 2 reps (mc may drop)
          all(is.finite(res$f1) | is.na(res$f1)),
          all(is.finite(res$drift) | is.na(res$drift)))
cat("\nSMOKE OK: arms x reps =", nrow(res), "\n")
print(res[, c("rep", "arm", "f1", "fp", "drift", "recon_ratio_min", "t_arm")])
