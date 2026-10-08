## 02e_outcome_na.R -- Outcome has 23 NA donors; test NA-handling modes with the raw (un-z-scored) reps.
here <- local({ f <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]); dirname(normalizePath(f)) })
source(file.path(here, "02b_metrics_variants.R"))
donors <- L$donors; meta <- sp$meta[donors, ]; reps <- make_reps(donors, "none", "none")
refv <- setNames(ref$value, paste(ref$representation, ref$covariate))
y0 <- as.character(meta$Outcome)
mode_f <- function(X, mode) {
  y <- y0
  if (mode == "drop") { ok <- !is.na(y); X <- X[ok, ]; y <- y[ok]; p <- class::knn.cv(X, factor(y), 5); return(macro_f1(y, p)) }
  y[is.na(y)] <- "NA_label"; p <- as.character(class::knn.cv(X, factor(y), 5))
  if (mode == "keep_all") return(macro_f1(y, p))                       # NA is a class, in the average
  ev <- y != "NA_label"                                                 # train on all, score non-NA donors only
  if (mode == "keep_train_eval_nonNA") return(macro_f1(y[ev], p[ev]))   # F1 over the 6 real classes
  if (mode == "keep_all_donors_F1_real_classes") { cl <- setdiff(unique(y), "NA_label"); return(mean(vapply(cl, function(c) f1_one(y, p, c), 0))) }
}
res <- sapply(c("drop", "keep_all", "keep_train_eval_nonNA", "keep_all_donors_F1_real_classes"), function(m)
  sapply(reps_scored, function(rn) mean(sapply(1:100, function(sd) { set.seed(sd); mode_f(reps[[rn]], m) }))))
res <- cbind(res, ref = refv[paste(reps_scored, "Outcome")])
print(round(res, 3)); cat("MAD vs ref:\n"); print(round(colMeans(abs(res[, 1:4] - res[, "ref"])), 4))
