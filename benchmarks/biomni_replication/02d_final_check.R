## 02d_final_check.R -- writes results/diag/knn_variants.csv (top 5 + headline + current) and tests the headline variant
## (raw, un-z-scored reps + class::knn.cv, i.e. the current pipeline minus all z-scoring) against the tie-noise floor.
here <- local({ f <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]); dirname(normalizePath(f)) })
source(file.path(here, "02b_metrics_variants.R"))
o <- readRDS(file.path(here, "results/diag/knn_variants.rds")); s <- subset(o$summ, f1 == "truth")
s$label <- with(s, sprintf("comp=%s pca=%s dist=%s knn=%s NA=%s donors=%s", comp_scale, pca_scale, dist, knn, na, donors))
top <- head(s[order(s$mean_abs_diff, s$n_departures), ], 5); top$role <- "top5_by_mad"
head_row <- subset(s, comp_scale == "none" & pca_scale == "none" & dist == "euclid" & knn == "random_useall" & na == "drop" & donors == "122")
head_row$role <- "headline_simplest(=class::knn.cv default, no z-scoring)"
cur <- subset(s, comp_scale == "z" & pca_scale == "zbefore" & dist == "euclid" & knn == "random_useall" & na == "drop" & donors == "122")
cur$role <- "current_02_metrics.R"
keep <- c("role", "label", "mean_abs_diff", "n_cells_lt_0.03", "n_departures")
write.csv(rbind(top, head_row, cur)[, keep], file.path(here, "results/diag/knn_variants.csv"), row.names = FALSE)
print(rbind(top, head_row, cur)[, keep], row.names = FALSE)

## headline variant via class::knn.cv over many seeds
donors <- L$donors; meta <- sp$meta[donors, ]; reps <- make_reps(donors, "none", "none")
refv <- setNames(ref$value, paste(ref$representation, ref$covariate))
cells <- expand.grid(representation = reps_scored, covariate = cats, stringsAsFactors = FALSE)
NS <- 200
sims <- sapply(seq_len(NS), function(sd) { set.seed(sd); mapply(function(rn, cv) knn_macro_f1(reps[[rn]], meta[[cv]], 5L), cells$representation, cells$covariate) })
cells$ref <- refv[paste(cells$representation, cells$covariate)]
cells$expected <- rowMeans(sims); cells$sd <- apply(sims, 1, sd); cells$seed1 <- sims[, 1]
cells$z <- (cells$ref - cells$expected) / pmax(cells$sd, 1e-9)
cells$absdiff_expected <- abs(cells$expected - cells$ref)
cells$pctl <- vapply(seq_len(nrow(cells)), function(i) mean(sims[i, ] <= cells$ref[i]), 0)
write.csv(cells, file.path(here, "results/diag/headline_cells.csv"), row.names = FALSE)
options(width = 200)
print(cells[, c("representation", "covariate", "ref", "expected", "sd", "seed1", "z")], digits = 3, row.names = FALSE)
cat("\nMAD of expectation vs ref, by rep:\n"); print(tapply(cells$absdiff_expected, cells$representation, mean), digits = 3)
cat("overall MAD (expected):", mean(cells$absdiff_expected), " seed1 draw:", mean(abs(cells$seed1 - cells$ref)), "\n")
## null: if headline were the true generator, MAD of a single draw vs the expectation
null <- apply(sims, 2, function(x) mean(abs(x - cells$expected)))
cat("tie-noise floor: MAD(single draw vs expectation) mean", mean(null), " 95% range", quantile(null, c(.025, .975)), "\n")
comp <- cells$representation == "composition"
cat("composition-only MAD (expected)", mean(cells$absdiff_expected[comp]), " floor", mean(apply(sims[comp, ], 2, function(x) mean(abs(x - cells$expected[comp])))), "\n")
