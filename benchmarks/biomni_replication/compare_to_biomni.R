## compare_to_biomni.R <our_eval_csv> [ref_eval_csv]  -- per-cell |diff| vs Biomni eval_all.csv, flags > 0.03.
a <- commandArgs(TRUE)
here <- local({ f <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]); dirname(normalizePath(f)) })
ours <- read.csv(a[1])
ref <- read.csv(if (length(a) > 1) a[2] else file.path(here, "ref/cellprograms_benchmark/combat/tables/eval_all.csv"))
if ("metric_variant" %in% names(ours)) ours <- ours[ours$metric_variant %in% c("biomni_zerofill", "biomni_spearman"), ]
m <- merge(ours[, c("representation", "covariate", "value")], ref[, c("representation", "covariate", "value")],
           by = c("representation", "covariate"), suffixes = c("_ours", "_biomni"))
m$diff <- m$value_ours - m$value_biomni; m$flag <- ifelse(abs(m$diff) > 0.03, "***", "")
m <- m[order(m$representation, m$covariate), ]
options(width = 160)
for (rn in unique(m$representation)) { cat("\n==", rn, "\n"); print(m[m$representation == rn, -1], row.names = FALSE, digits = 3) }
cat(sprintf("\n%d matched cells; %d flagged (|diff| > 0.03); mean |diff| = %.3f\n", nrow(m), sum(m$flag != ""),
            mean(abs(m$diff))))
cat("not in ours:", paste(setdiff(unique(ref$representation), ours$representation), collapse = ", "), "\n")
