## 09_effect_size.R : recompute the COMBAT sharing spectra (with z_pair / excess_frac) from the saved fits of 07
library(stats); library(utils); for (f in list.files("../../R", full.names = TRUE)) source(f)
for (n in c("none", "Institute", "Institute_Outcome")) {
  fit <- readRDS(sprintf("results/sharing/combat_repo_%s_fit.rds", n))
  ss <- suppressWarnings(sharing_spectrum(fit, space = "scores", n_perm = 1999, seed = 1))$summary
  write.csv(ss, sprintf("results/sharing/combat_repo_%s_pairs_v2.csv", n), row.names = FALSE)
}
