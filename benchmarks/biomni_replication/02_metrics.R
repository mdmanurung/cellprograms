## 02_metrics.R -- Biomni-style evaluation (reconstructed; their metrics.R/baselines.R were unavailable).
## Library: source("02_metrics.R").  CLI: Rscript 02_metrics.R <cohort> <fit_rds|label=rds[,..]|NONE> <out_csv>
## Metrics (Biomni): categorical = leave-one-out 5-NN macro-F1 (Euclidean, class::knn.cv), macro over
## classes PRESENT in the truth; continuous = max |Spearman| over representation columns.
## Missing donors are zero-filled AFTER per-column z-scoring (align_donors).

.script_dir <- function() {
  f <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
  if (length(f)) dirname(normalizePath(f[1])) else getwd()
}
DATA_ROOT <- Sys.getenv("BIOMNI_DATA", file.path(.script_dir(), "data"))

## ---------------------------------------------------------------- primitives
.zscore <- function(M) {
  M <- as.matrix(M)
  mu <- colMeans(M, na.rm = TRUE)
  s <- apply(M, 2, stats::sd, na.rm = TRUE); s[!is.finite(s) | s < 1e-12] <- 1
  M <- sweep(sweep(M, 2, mu, "-"), 2, s, "/")
  M[is.na(M)] <- 0
  M
}

align_donors <- function(sc, donors) {
  out <- matrix(0, length(donors), ncol(sc), dimnames = list(donors, colnames(sc)))
  common <- intersect(rownames(sc), donors)
  out[common, ] <- sc[common, , drop = FALSE]
  out
}

## PCA scores (k comps, NA -> column mean); raw (unscaled) scores, z-scoring happens in blocked_rep
.safe_pca <- function(X, k) {
  X <- as.matrix(X)
  X <- X[, apply(X, 2, function(v) sum(is.finite(v)) > 0), drop = FALSE]
  cm <- colMeans(X, na.rm = TRUE)
  for (j in which(colSums(is.na(X)) > 0)) X[is.na(X[, j]), j] <- cm[j]
  X <- X[, apply(X, 2, stats::sd) > 1e-12, drop = FALSE]
  k <- min(k, ncol(X), nrow(X) - 1L)
  s <- stats::prcomp(X, center = TRUE, scale. = FALSE, rank. = k)
  out <- s$x[, seq_len(k), drop = FALSE]
  colnames(out) <- paste0("PC", seq_len(k)); out
}

rep_random <- function(donors, k = 10L, seed = 1L) {
  set.seed(seed)
  matrix(stats::rnorm(length(donors) * k), length(donors), k, dimnames = list(donors, paste0("R", seq_len(k))))
}
## Baseline scaling: Biomni-faithful = NO z-scoring (reproduces its COMBAT baselines: global/perct/grouped PCA Death28/sex/Institute
## match to 3 decimals, continuous max|Spearman| exactly). Set options(cp.baseline_zscore = TRUE) for the old behaviour.
.bz <- function(M) if (isTRUE(getOption("cp.baseline_zscore", FALSE))) .zscore(M) else as.matrix(M)
rep_composition <- function(comp) { p <- comp / rowSums(comp); .bz(p) }
rep_clr <- function(comp) {
  p <- comp / rowSums(comp)
  ## counts -> +0.5 pseudocount; proportions -> +1e-4 (ponytail: Biomni's pseudocount unknown)
  x <- if (all(abs(comp - round(comp)) < 1e-8)) comp + 0.5 else p + 1e-4
  l <- log(x / rowSums(x)); .bz(l - rowMeans(l))
}

## Blocked representation: list of block name -> (observed donors x cols) matrices.
## z-scores each column on observed donors, zero-fills the rest, remembers block structure.
blocked_rep <- function(lst, donors, zscore = TRUE) {
  lst <- lst[vapply(lst, function(m) !is.null(m) && ncol(m) > 0L, logical(1))]
  mats <- lapply(names(lst), function(b) {
    m <- lst[[b]]; m <- m[rownames(m) %in% donors, , drop = FALSE]
    align_donors(if (zscore) .zscore(m) else m, donors)
  })
  X <- do.call(cbind, mats)
  blocks <- rep(names(lst), vapply(mats, ncol, integer(1)))
  if (anyDuplicated(colnames(X))) colnames(X) <- paste(blocks, colnames(X), sep = "__")
  obs <- vapply(lst, function(m) donors %in% rownames(m), logical(length(donors)))
  dimnames(obs) <- list(donors, names(lst))
  attr(X, "blocks") <- blocks; attr(X, "observed") <- obs
  X
}
single_block <- function(X, name = "all") {  # fully observed, one block
  attr(X, "blocks") <- rep(name, ncol(X))
  attr(X, "observed") <- matrix(TRUE, nrow(X), 1L, dimnames = list(rownames(X), name)); X
}

## ---------------------------------------------------------------- cohort I/O
AGE_MAP <- c("19-30" = 24.5, "31-40" = 35.5, "41-50" = 45.5, "51-60" = 55.5,
             "61-70" = 65.5, "71-80" = 75.5, "81-90" = 85.5, ">=91" = 92.5)

.find <- function(d, f) { p <- file.path(c(d, file.path(d, "mat")), f); p[file.exists(p)][1] }
.rd <- function(f) as.matrix(utils::read.csv(f, row.names = 1, check.names = FALSE))

load_cohort <- function(cohort, root = DATA_ROOT) {
  d <- file.path(root, cohort)
  cf <- Sys.glob(file.path(d, "mat", "ct_*.csv"))
  mats <- lapply(cf, .rd)
  names(mats) <- sub("\\.csv$", "", sub("^ct_", "", basename(cf)))
  gf <- c(Sys.glob(file.path(d, "group_*.csv")), Sys.glob(file.path(d, "mat", "group_*.csv")))
  groups <- lapply(gf, .rd); names(groups) <- sub("\\.csv$", "", sub("^group_", "", basename(gf)))
  meta <- utils::read.csv(file.path(d, "donor_meta.csv"), row.names = 1, check.names = FALSE,
                          na.strings = c("NA", ""))
  list(cohort = cohort, mats = mats, groups = groups,
       global = .rd(.find(d, "global.csv")), comp = .rd(.find(d, "composition.csv")),
       meta = meta, donors = sort(unique(unlist(lapply(mats, rownames)))))
}

## Covariate spec. COMBAT: exactly Biomni's. Other cohorts: heuristic. Candidates are character/factor columns
## with 2-12 levels plus numeric columns with >= 2 distinct values.
## TECHNICAL covariates (Panel B; lower = better, never counted as biology) are recognised by NAME, for
## categorical AND continuous columns (OneK1K pool_number is numeric):
##   pool_number/pool, Institute, Site, GEX_region/region, sequencing_platform/platform, study, batch,
##   chemistry, assay, lane, Collection_Day, Resample, dataset, tissue_source (see TECH_REGEX).
## Numeric BIOLOGY covariates additionally need >= 6 distinct values. The per-cohort biology/technical/
## excluded lists used for the verdict tally live in verdicts_cohort_spec.R (summarize_eval.R).
TECH_REGEX <- paste0("batch|site|institute|platform|region|study|chemistry|assay|lane|pool|",
                     "collection_?day|resample|dataset|tissue_?source")

cohort_spec <- function(meta, cohort = "combat") {
  if (grepl("combat", cohort, ignore.case = TRUE)) {
    if (is.character(meta$Age) || is.factor(meta$Age)) meta$Age_num <- unname(AGE_MAP[as.character(meta$Age)])
    else meta$Age_num <- as.numeric(meta$Age)
    d28 <- as.character(meta$Death28)
    meta$Death28 <- ifelse(is.na(d28), NA, d28 %in% c("True", "TRUE", "true", "1"))  # Biomni: == "True" (NA kept)
    return(list(meta = meta, cat_bio = c("Source", "Outcome", "Death28", "sex"),
                cont_bio = c("Age_num", "Hospitalstay", "TimeSinceOnset"),
                cat_tech = c("Institute", "GEX_region"), cont_tech = character(0)))
  }
  isnum <- vapply(meta, is.numeric, logical(1))
  nl <- vapply(meta, function(v) length(unique(stats::na.omit(v))), integer(1))
  is_tech <- grepl(TECH_REGEX, names(meta), ignore.case = TRUE)
  cat_all <- !isnum & nl >= 2 & nl <= 12
  list(meta = meta,
       cat_bio = names(meta)[cat_all & !is_tech], cat_tech = names(meta)[cat_all & is_tech],
       cont_bio = names(meta)[isnum & nl >= 6 & !is_tech], cont_tech = names(meta)[isnum & nl >= 2 & is_tech])
}

## ---------------------------------------------------------------- representations
build_baselines <- function(L, seed = 1L, pca_zscore = isTRUE(getOption("cp.baseline_zscore", FALSE))) {
  donors <- L$donors
  reps <- list()
  reps$random <- single_block(rep_random(donors, 10L, seed), "random")
  reps$composition <- single_block(align_donors(rep_composition(L$comp), donors), "composition")
  reps$clr_composition <- single_block(align_donors(rep_clr(L$comp), donors), "clr_composition")
  reps$global_pca <- blocked_rep(list(global = .safe_pca(L$global, 10L)), donors, pca_zscore)
  reps$perct_pca <- blocked_rep(lapply(L$mats, .safe_pca, k = 3L), donors, pca_zscore)
  if (length(L$groups)) reps$grouped_pca <- blocked_rep(lapply(L$groups, .safe_pca, k = 5L), donors, pca_zscore)  # grouped matrices only built for COMBAT
  reps
}

## cp representation from a fitted + canonicalized fit: per-CT program scores, z-scored over observed
## donors, zero-filled for donors lacking that CT.
cp_representation <- function(fit, donors) {
  sc <- if (requireNamespace("cellprograms", quietly = TRUE)) cellprograms::program_scores(fit, "list") else fit$scores
  blocked_rep(sc, donors, zscore = TRUE)
}

parse_fits <- function(arg) {  # "NONE" | "path" | "label=path,label2=path2"
  if (is.null(arg) || toupper(arg) == "NONE") return(list())
  parts <- strsplit(arg, ",")[[1]]
  lab <- ifelse(grepl("=", parts), sub("=.*", "", parts), "cp_L0")
  fits <- lapply(sub(".*=", "", parts), function(p) { f <- readRDS(p); f })
  stats::setNames(fits, lab)
}

all_representations <- function(L, fits, ...) {
  reps <- build_baselines(L, ...)
  for (nm in names(fits)) reps[[nm]] <- cp_representation(fits[[nm]], L$donors)
  reps
}

## ---------------------------------------------------------------- metrics
macro_f1 <- function(truth, pred) {  # average over classes present in the truth (fixed Biomni bug)
  truth <- as.character(truth); pred <- as.character(pred)
  mean(vapply(unique(truth), function(cl) {
    d <- sum(truth == cl) + sum(pred == cl)
    if (d == 0) 0 else 2 * sum(truth == cl & pred == cl) / d
  }, numeric(1)))
}

knn_macro_f1 <- function(X, y, k = 5L) {
  ok <- !is.na(y)
  if (sum(ok) < k + 2 || length(unique(y[ok])) < 2) return(NA_real_)
  pred <- class::knn.cv(unclass_matrix(X)[ok, , drop = FALSE], factor(y[ok]), k = k)
  macro_f1(y[ok], pred)
}
unclass_matrix <- function(X) { attributes(X) <- attributes(X)[c("dim", "dimnames")]; X }

max_abs_spearman <- function(X, y) {
  ok <- !is.na(y)
  if (sum(ok) < 5 || stats::sd(y[ok]) == 0) return(NA_real_)
  r <- suppressWarnings(stats::cor(unclass_matrix(X)[ok, , drop = FALSE], y[ok], method = "spearman"))
  if (all(is.na(r))) NA_real_ else max(abs(r), na.rm = TRUE)
}

evaluate_representation <- function(rep, meta, cat_cols, cont_cols, k = 5L, seed = 1L) {
  set.seed(seed)  # class::knn breaks ties at random
  cat_cols <- intersect(cat_cols, names(meta)); cont_cols <- intersect(cont_cols, names(meta))
  rows <- c(
    lapply(cat_cols, function(cv) data.frame(covariate = cv, type = "categorical", metric = "knn_macro_f1",
                                             value = knn_macro_f1(rep, meta[[cv]], k))),
    lapply(cont_cols, function(cv) data.frame(covariate = cv, type = "continuous", metric = "max_abs_spearman",
                                              value = max_abs_spearman(rep, as.numeric(meta[[cv]])))))
  do.call(rbind, rows)
}

eval_all <- function(reps, spec, donors) {
  meta_al <- spec$meta[donors, , drop = FALSE]  # donors absent from meta -> all-NA rows
  rownames(meta_al) <- donors
  out <- do.call(rbind, lapply(names(reps), function(rn) {
    ev <- evaluate_representation(reps[[rn]], meta_al, c(spec$cat_bio, spec$cat_tech),
                                  c(spec$cont_bio, spec$cont_tech))
    data.frame(representation = rn, ev)
  }))
  out$panel <- ifelse(out$covariate %in% c(spec$cat_tech, spec$cont_tech), "B_technical", "A_biology")
  rownames(out) <- NULL
  out
}

## ---------------------------------------------------------------- CLI
if (!interactive() && identical(basename(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])),
                                "02_metrics.R")) {
  a <- commandArgs(TRUE)
  if (length(a) != 3) stop("usage: Rscript 02_metrics.R <cohort> <fit_rds|NONE> <out_csv>")
  fits <- parse_fits(a[2])
  if (length(fits)) pkgload::load_all(file.path(.script_dir(), "fork"), quiet = TRUE)
  L <- load_cohort(a[1]); spec <- cohort_spec(L$meta, a[1])
  res <- eval_all(all_representations(L, fits), spec, L$donors)
  utils::write.csv(res, a[3], row.names = FALSE)
  cat("wrote", nrow(res), "rows to", a[3], "\n")
}
