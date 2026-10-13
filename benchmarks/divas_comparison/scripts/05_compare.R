# 05_compare.R -- head-to-head comparison: DIVAS vs cellprograms
#
# Inputs (out_dir of scripts 02-04):
#   divas6_fit.rds, cmp_marginal.rds, sharing_geneblocks.rds  (cellprograms)
#   divas_scores.csv, divas_components.csv, divas_top100_loadings.csv,
#   divas_loadings.rds                                        (DIVAS)
#   input/final_sample_info.csv                               (clinical metadata)
#
# Outputs: comparison tables + figures (PNG) in out_dir.
#
# Usage: Rscript 05_compare.R <data_dir> <out_dir>

suppressPackageStartupMessages({
  if (requireNamespace("cellprograms", quietly = TRUE)) library(cellprograms)
  else pkgload::load_all(Sys.getenv("CELLPROGRAMS_ROOT", unset = file.path("..", "..", "..")),
                         quiet = TRUE)
  library(ggplot2)
})

args <- commandArgs(trailingOnly = TRUE)
data_dir <- if (length(args) >= 1) args[1] else "data"
out_dir  <- if (length(args) >= 2) args[2] else "results"

BLOCK_LABEL <- c(CD4_T = "CD4", CD8_T = "CD8", CD14_Monocyte = "CD14", NK = "NK",
                 proteomics = "pro", metabolomics = "metab")

fit   <- readRDS(file.path(out_dir, "divas6_fit.rds"))
cmp   <- readRDS(file.path(out_dir, "cmp_marginal.rds"))
share <- readRDS(file.path(out_dir, "sharing_geneblocks.rds"))

dv_scores <- read.csv(file.path(out_dir, "divas_scores.csv"), row.names = 1,
                      check.names = FALSE)
dv_annot  <- read.csv(file.path(out_dir, "divas_components.csv"),
                      stringsAsFactors = FALSE)
dv_top    <- read.csv(file.path(out_dir, "divas_top100_loadings.csv"),
                      stringsAsFactors = FALSE)

meta <- read.csv(file.path(data_dir, "input", "final_sample_info.csv"),
                 stringsAsFactors = FALSE)
meta <- meta[match(rownames(dv_scores), meta$sample_id), ]
stopifnot(!anyNA(meta$sample_id))
wos <- meta$wos
sev <- factor(meta$severity, levels = c("mild", "moderate", "severe"))

# ---- cellprograms score matrix: samples x programs (short labels: CD4_1, ...)
cp_scores <- do.call(cbind, lapply(fit$cell_types, function(ct) {
  S <- fit$scores[[ct]]
  colnames(S) <- paste0(BLOCK_LABEL[[ct]], "_",
                        sub(paste0("^", ct, "_"), "", colnames(S)))
  S
}))
stopifnot(identical(rownames(cp_scores), rownames(dv_scores)))

# ---- CMP consensus scores, if the cmp result carries them
cmp_scores <- NULL
for (nm in c("consensus_scores", "scores")) {
  cand <- cmp[[nm]]
  if (!is.null(dim(cand)) && nrow(cand) == nrow(dv_scores) &&
      !identical(colnames(cand), colnames(cp_scores))) {
    cmp_scores <- cand
    break
  }
}
cp_mat <- if (!is.null(cmp_scores)) cbind(cp_scores, cmp_scores) else cp_scores
cp_mat <- cp_mat[, make.unique(colnames(cp_mat)), drop = FALSE]
# constant score columns carry no concordance information and poison cor()
const_cp <- apply(cp_mat, 2, sd) == 0
if (any(const_cp)) {
  warning(sprintf("dropping %d constant cellprograms score column(s)", sum(const_cp)))
  cp_mat <- cp_mat[, !const_cp, drop = FALSE]
}
cat("cellprograms scores:", ncol(cp_scores), "programs +",
    if (!is.null(cmp_scores)) ncol(cmp_scores) else 0, "CMPs\n")

spcor <- function(x, y) {
  suppressWarnings(cor(x, y, method = "spearman", use = "pairwise.complete.obs"))
}

# ================= 1. score concordance =================
cat("== score concordance\n")
cor_mat <- matrix(NA_real_, ncol(dv_scores), ncol(cp_mat),
                  dimnames = list(colnames(dv_scores), colnames(cp_mat)))
for (j in seq_len(ncol(dv_scores))) {
  for (k in seq_len(ncol(cp_mat))) {
    cor_mat[j, k] <- abs(spcor(dv_scores[, j], cp_mat[, k]))
  }
}
# NA cor (constant inputs) means no concordance evidence -> score as 0
n_na <- sum(is.na(cor_mat))
if (n_na) warning(sprintf("%d NA correlations scored as 0", n_na))
cor_mat[is.na(cor_mat)] <- 0
cmp_cols <- if (!is.null(cmp_scores)) colnames(cmp_scores) else character(0)
concord <- data.frame(
  divas_component = colnames(dv_scores),
  divas_class = dv_annot$divas_class,
  divas_pattern = dv_annot$pattern,
  best_cp_program = colnames(cp_mat)[max.col(cor_mat)],
  best_abs_cor = apply(cor_mat, 1, max),
  n_cp_above_0.5 = rowSums(cor_mat > 0.5),
  best_is_cmp = colnames(cp_mat)[max.col(cor_mat)] %in% cmp_cols,
  stringsAsFactors = FALSE
)
write.csv(concord, file.path(out_dir, "concordance_divas_to_cellprograms.csv"),
          row.names = FALSE)
write.csv(round(cor_mat, 4), file.path(out_dir, "score_cor_matrix.csv"))

rev <- data.frame(
  cp_program = colnames(cp_mat),
  best_divas_component = colnames(dv_scores)[max.col(t(cor_mat))],
  best_abs_cor = apply(cor_mat, 2, max),
  stringsAsFactors = FALSE
)
write.csv(rev, file.path(out_dir, "concordance_cellprograms_to_divas.csv"),
          row.names = FALSE)

cat(sprintf("DIVAS components with a cellprograms counterpart (|cor|>0.5): %d/%d\n",
            sum(concord$best_abs_cor > 0.5), nrow(concord)))
cat(sprintf("cellprograms scores with a DIVAS counterpart (|cor|>0.5): %d/%d\n",
            sum(rev$best_abs_cor > 0.5), nrow(rev)))

# ================= 2. loading overlap for mutual-best gene-block pairs =================
cat("== loading overlap\n")
cp_top <- do.call(rbind, lapply(fit$cell_types, function(ct) {
  W <- fit$loadings[[ct]]
  lab <- BLOCK_LABEL[[ct]]
  do.call(rbind, lapply(seq_len(ncol(W)), function(k) {
    o <- order(-abs(W[, k]))[seq_len(min(100L, nrow(W)))]
    data.frame(block = lab, program = paste0(lab, "_", k), rank = seq_along(o),
               feature = rownames(W)[o], stringsAsFactors = FALSE)
  }))
}))
# mutual best matches: forward and reverse argmax agree
fwd_prog <- setNames(concord$best_cp_program, concord$divas_component)
rev_comp <- setNames(rev$best_divas_component, rev$cp_program)
mutual <- names(fwd_prog)[names(fwd_prog) == unname(rev_comp[unname(fwd_prog)])]  # comp -> prog -> comp round-trip
mutual <- mutual[!is.na(mutual)]

ov <- lapply(mutual, function(comp) {
  prog <- unname(fwd_prog[[comp]])
  prog_lab <- sub("_[0-9]+$", "", prog)
  dv_blk <- dv_top[dv_top$component == comp & dv_top$block == prog_lab, ]
  cp_blk <- cp_top[cp_top$program == prog, ]
  if (nrow(dv_blk) == 0L || nrow(cp_blk) == 0L) return(NULL)
  data.frame(divas_component = comp, cp_program = prog, block = prog_lab,
             jaccard_top100 =
               length(intersect(dv_blk$feature, cp_blk$feature)) /
               length(union(dv_blk$feature, cp_blk$feature)),
             stringsAsFactors = FALSE)
})
overlap <- do.call(rbind, ov)
write.csv(overlap, file.path(out_dir, "loading_overlap_top100.csv"), row.names = FALSE)
cat(sprintf("top-100 loading Jaccard for %d mutual-best gene-block pairs: median %.2f (IQR %.2f-%.2f)\n",
            nrow(overlap), median(overlap$jaccard_top100),
            quantile(overlap$jaccard_top100, 0.25), quantile(overlap$jaccard_top100, 0.75)))

# ================= 3. taxonomy agreement =================
cat("== taxonomy agreement\n")
cp_class <- share$summary[, c("program_id", "cell_type", "classification", "n_conditions")]
cp_class$cp_program <- paste0(BLOCK_LABEL[cp_class$cell_type], "_",
                              sub(".*_", "", cp_class$program_id))  # strip "<celltype>_" prefix
tax <- merge(cp_class[, c("cp_program", "classification", "n_conditions")], rev,
             by = "cp_program")
tax <- merge(tax, dv_annot[, c("component", "n_blocks", "divas_class")],
             by.x = "best_divas_component", by.y = "component", all.x = TRUE)
tax$divas_breadth <- ifelse(is.na(tax$divas_class), NA,
                     ifelse(tax$n_blocks == 6, "joint",
                     ifelse(tax$n_blocks == 1, "individual", "partial")))
tax$cp_breadth <- c(shared = "shared", private = "private", partial = "partial")[tax$classification]
tax_tab <- table(cellprograms = factor(tax$cp_breadth,
                                       levels = c("shared", "partial", "private")),
                 divas = factor(tax$divas_breadth,
                                levels = c("joint", "partial", "individual")))
cat("agreement cross-tab (matched gene-block programs):\n"); print(tax_tab)
write.csv(as.data.frame(tax_tab), file.path(out_dir, "taxonomy_agreement.csv"),
          row.names = FALSE)

# ================= 4. clinical association =================
cat("== clinical association\n")
assoc <- function(S, label) {
  pw <- apply(S, 2, function(v) suppressWarnings(cor(v, wos, method = "spearman",
                                                     use = "pairwise.complete.obs")))
  pp <- apply(S, 2, function(v)
    suppressWarnings(cor.test(v, wos, method = "spearman",
                              exact = FALSE)$p.value))
  kw <- apply(S, 2, function(v) suppressWarnings(kruskal.test(v ~ sev)$p.value))
  data.frame(method = label, score = colnames(S), rho_wos = pw, p_wos = pp,
             p_severity = kw, stringsAsFactors = FALSE)
}
a_div <- assoc(dv_scores, "DIVAS")
a_cp  <- assoc(cp_mat, "cellprograms")
a_all <- rbind(a_div, a_cp)
a_all$p_wos_bh <- unlist(lapply(c("DIVAS", "cellprograms"), function(m)
  pmin(1, p.adjust(a_all$p_wos[a_all$method == m], method = "BH"))))
a_all$p_sev_bh <- unlist(lapply(c("DIVAS", "cellprograms"), function(m)
  pmin(1, p.adjust(a_all$p_severity[a_all$method == m], method = "BH"))))
write.csv(a_all, file.path(out_dir, "clinical_associations.csv"), row.names = FALSE)
cat("wos associations at BH<0.05: DIVAS",
    sum(a_all$p_wos_bh < 0.05 & a_all$method == "DIVAS"), "/",
    sum(a_all$method == "DIVAS"), "| cellprograms",
    sum(a_all$p_wos_bh < 0.05 & a_all$method == "cellprograms"), "/",
    sum(a_all$method == "cellprograms"), "\n")

# cross-validated ridge prediction of wos (identical folds for all methods)
cat("== CV ridge prediction of wos\n")
set.seed(1)
folds <- sample(rep(seq_len(10), length.out = nrow(dv_scores)))
cv_ridge <- function(S, y, folds) {
  S <- as.matrix(S)
  pred <- numeric(length(y))
  for (f in seq_len(10)) {
    tr <- folds != f
    te <- folds == f
    Xtr <- scale(S[tr, , drop = FALSE])
    mu <- attr(Xtr, "scaled:center")
    sd <- attr(Xtr, "scaled:scale")
    sd[!is.finite(sd) | sd == 0] <- 1
    Xte <- scale(S[te, , drop = FALSE], center = mu, scale = sd)
    fit_cv <- glmnet::cv.glmnet(Xtr, y[tr], alpha = 0, nfolds = 5)
    pred[te] <- as.numeric(predict(fit_cv, newx = Xte, s = "lambda.min"))
  }
  c(R2 = cor(y, pred)^2, RMSE = sqrt(mean((y - pred)^2)))
}
res_div  <- cv_ridge(dv_scores, wos, folds)
res_cp   <- cv_ridge(cp_mat, wos, folds)
sev_dummies <- model.matrix(~ sev)[, -1, drop = FALSE]  # moderate + severe indicators
res_base <- cv_ridge(sev_dummies, wos, folds)
cv_tab <- data.frame(
  predictors = c("DIVAS components", "cellprograms scores", "severity-only baseline"),
  n_predictors = c(ncol(dv_scores), ncol(cp_mat), ncol(sev_dummies)),
  R2 = c(res_div[["R2"]], res_cp[["R2"]], res_base[["R2"]]),
  RMSE = c(res_div[["RMSE"]], res_cp[["RMSE"]], res_base[["RMSE"]])
)
write.csv(cv_tab, file.path(out_dir, "cv_wos_prediction.csv"), row.names = FALSE)
print(cv_tab)

# ================= 5. unique findings =================
cat("== unique components/programs\n")
uniq_div <- concord[concord$best_abs_cor < 0.3, ]
uniq_cp  <- rev[rev$best_abs_cor < 0.3, ]
write.csv(uniq_div, file.path(out_dir, "unique_divas_components.csv"), row.names = FALSE)
write.csv(uniq_cp, file.path(out_dir, "unique_cellprograms_programs.csv"), row.names = FALSE)
cat("unique (no counterpart, |cor|<0.3): DIVAS", nrow(uniq_div), "components; cellprograms",
    nrow(uniq_cp), "scores\n")

# ================= figures =================
cat("== figures\n")
pal <- c("#0072B2", "#E69F00", "#009E73", "#CC79A7", "#56B4E9", "#D55E00")
theme_set(theme_minimal(base_family = "Liberation Sans"))

p1 <- ggplot(concord, aes(x = best_abs_cor, fill = divas_class)) +
  geom_histogram(binwidth = 0.05, boundary = 0, color = "white", linewidth = 0.1) +
  geom_vline(xintercept = 0.5, linetype = "dashed", color = "grey30") +
  scale_fill_manual(values = c(joint = pal[3], partial = pal[2], individual = pal[1]),
                    name = "DIVAS class") +
  labs(x = "Best |Spearman cor| to a cellprograms score",
       y = "DIVAS components",
       title = "Score concordance: DIVAS vs cellprograms")
ggsave(file.path(out_dir, "fig_concordance.png"), p1, width = 7, height = 4.2, dpi = 300)

tax_df <- as.data.frame(tax_tab)
if (sum(tax_df$Freq) > 0) {
  p2 <- ggplot(tax_df, aes(x = divas, y = cellprograms, fill = Freq)) +
    geom_tile(color = "white") +
    geom_text(aes(label = Freq), color = "grey15", size = 4) +
    scale_fill_gradient(low = "#F5F5F5", high = pal[1], name = "n programs") +
    labs(x = "DIVAS pattern breadth (matched component)",
         y = "cellprograms sharing class",
         title = "Taxonomy agreement, matched gene-block programs")
  ggsave(file.path(out_dir, "fig_taxonomy.png"), p2, width = 6, height = 4, dpi = 300)
}

p3 <- ggplot(a_all, aes(x = rho_wos, y = -log10(pmax(p_wos_bh, 1e-300)), color = method)) +
  geom_point(alpha = 0.75, size = 1.6) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey40") +
  scale_color_manual(values = c(DIVAS = pal[1], cellprograms = pal[2])) +
  labs(x = "Spearman rho with WHO ordinal score",
       y = "-log10(BH-adjusted p)", color = "Method",
       title = "Clinical association of components/programs")
ggsave(file.path(out_dir, "fig_clinical.png"), p3, width = 7, height = 4.2, dpi = 300)

p4 <- ggplot(cv_tab, aes(x = reorder(predictors, R2), y = R2, fill = predictors)) +
  geom_col(width = 0.7) +
  geom_text(aes(label = sprintf("%.2f", R2)), hjust = -0.15, size = 3.6) +
  scale_fill_manual(values = pal[c(1, 2, 4)]) +
  coord_flip(ylim = c(0, max(cv_tab$R2) * 1.2)) +
  labs(x = NULL, y = "Cross-validated R2 (10-fold, ridge)",
       title = "Prediction of WHO ordinal score") +
  theme(legend.position = "none")
ggsave(file.path(out_dir, "fig_cv_prediction.png"), p4, width = 6.5, height = 3.6, dpi = 300)

cat("DONE\n")
