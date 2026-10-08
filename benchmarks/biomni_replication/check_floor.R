## check_floor.R -- null distribution of max|Spearman| at n = 122 donors (independent Gaussian X and y).
## Expected mean: ~0.17 at 10 dims, ~0.21 at 30, ~0.27 at 266 (tolerance 0.02). Base R only.
set.seed(42); n <- 122L; nrep <- 300L
target <- c("10" = 0.17, "30" = 0.21, "266" = 0.27)
ok <- TRUE
for (p in names(target)) {
  m <- replicate(nrep, max(abs(cor(matrix(rnorm(n * as.integer(p)), n), rnorm(n), method = "spearman"))))
  pass <- abs(mean(m) - target[[p]]) <= 0.02
  ok <- ok && pass
  cat(sprintf("p=%3s  mean max|rho| = %.3f (sd %.3f)  target %.2f  %s\n", p, mean(m), sd(m), target[[p]],
              if (pass) "PASS" else "FAIL"))
}
if (!ok) quit(status = 1)
