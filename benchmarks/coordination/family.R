## family.R -- name the comparison family of the primary decision table.
## prereg-coord-v2 reports the arm-versus-baseline contrasts WITHOUT multiplicity adjustment
## (user decision 2026-11-17), so the table must say which contrasts form the family and how
## many there are, and give the unadjusted two-sided t p-value of each, so a reader can adjust.
add_family <- function(paired, label) {
  if (!nrow(paired)) return(paired)
  se <- (paired$ci_hi - paired$mean_diff) / 1.96        # ci_hi = mean + 1.96 * sd / sqrt(n)
  tstat <- ifelse(se > 0, paired$mean_diff / se, NA_real_)
  paired$p_unadjusted <- 2 * stats::pt(-abs(tstat), df = paired$n_pairs - 1)
  paired$family <- label
  paired$family_size <- nrow(paired)
  paired$adjustment <- "none"
  paired
}
