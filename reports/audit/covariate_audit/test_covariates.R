#!/usr/bin/env Rscript
# test_covariates.R -- adversarial audit of cellprograms covariate functionality
# Run: R_LIBS=/workspace/.Rlib /opt/conda/bin/Rscript test_covariates.R
suppressPackageStartupMessages(library(cellprograms))

cat("== cellprograms covariate audit ==\n")
cat("cellprograms ", as.character(packageVersion("cellprograms")),
    "| flashier ", as.character(packageVersion("flashier")),
    "| flashier.utils ", as.character(packageVersion("flashier.utils")), "\n")
cat(R.version.string, "\n\n")

PASS <- 0L; FAIL <- 0L
ok <- function(id, cond, detail = "") {
  if (isTRUE(cond)) { PASS <<- PASS + 1L; cat(sprintf("[PASS] %s  %s\n", id, detail)) }
  else              { FAIL <<- FAIL + 1L; cat(sprintf("[FAIL] %s  %s\n", id, detail)) }
  invisible(NULL)
}
# capture value + error + warnings + messages in one pass
grab <- function(expr) {
  warns <- character(0); msgs <- character(0); err <- NULL; val <- NULL
  val <- tryCatch(
    withCallingHandlers(
      expr,
      warning = function(w) { warns <<- c(warns, conditionMessage(w)); invokeRestart("muffleWarning") },
      message = function(m) { msgs <<- c(msgs, conditionMessage(m)); invokeRestart("muffleMessage") }
    ),
    error = function(e) { err <<- conditionMessage(e); NULL }
  )
  list(value = val, error = err, warnings = warns, messages = msgs)
}
cd <- cellprograms:::.cov_design

## ---------------- data ----------------
set.seed(7)
sim <- simulate_pseudobulk(scenario = "mixed", n_samples = 40, n_genes = 150,
                           program_size = 20, seed = 11)
pb <- sim$pseudobulk
ids <- rownames(pb$B)
age   <- rnorm(40)
sex   <- factor(rep(c("F", "M"), each = 20))
batch <- factor(rep(c("b1", "b2", "b3", "b4"), times = 10))
site  <- factor(c(rep("S1", 30), rep("S2", 10)))   # constant within obs1-30
meta <- data.frame(observation_id = ids, age = age, sex = sex, batch = batch,
                   site = site, age_copy = age, site_const = factor(rep("S1", 40)))
# plant a strong age-driven module in B (20 genes), to test removal
w_age <- rnorm(20)
pb$B[, 1:20] <- pb$B[, 1:20] + 3 * (age %*% matrix(w_age, nrow = 1))
x0 <- as_cell_program_data(pb, sample_metadata = meta)

FA <- list(features = "variable", n_variable_genes = 100, max_factors = 3,
           backfit = FALSE, feature_scale = "logcounts")
fit_none   <- do.call(fit_celltype_programs, c(list(x = x0), FA))
fit_res    <- do.call(fit_celltype_programs, c(list(x = x0, covariates = ~ age + sex + batch), FA))
fit_res_c  <- do.call(fit_celltype_programs, c(list(x = x0, covariates = c("age","sex","batch")), FA))
fit_fix    <- do.call(fit_celltype_programs, c(list(x = x0, covariates = ~ age + sex + batch,
                                                    covariate_mode = "fixed"), FA))
fit_collin <- do.call(fit_celltype_programs, c(list(x = x0, covariates = c("age","age_copy")), FA))

## ---------------- T1: ordering (selection before adjustment; residualize before center/scale)
sel_raw <- cellprograms:::.select_features(x0, "variable", 100, NULL)$B
ok("T1a variable-gene selection uses RAW Y (before adjustment)",
   identical(fit_res$preprocessing$feature_sets$B, sel_raw),
   sprintf("identical panels: %s", identical(fit_res$preprocessing$feature_sets$B, sel_raw)))
Yraw <- as.matrix(pb$B)[, fit_res$preprocessing$feature_sets$B, drop = FALSE]
Xres <- cd(rownames(Yraw), meta, ~ age + sex + batch, "B")
R    <- lm.fit(Xres, Yraw)$residuals
rp   <- fit_res$fits$B$preproc_replay
d1 <- max(abs(sweep(R, 2, rp$center, "-") - fit_res$fits$B$Y_used))
ok("T1b Y_used == center(residualize(Y)) i.e. residualize BEFORE centering", d1 < 1e-8,
   sprintf("max|diff| = %.2e", d1))
fit_res_s <- do.call(fit_celltype_programs,
                     c(list(x = x0, covariates = ~ age + sex + batch, scale = TRUE), FA))
rp_s <- fit_res_s$fits$B$preproc_replay
sd_res <- apply(R, 2, sd); sd_raw <- apply(Yraw, 2, sd)
ok("T1c with scale=TRUE, scale values are RESIDUAL sds (not raw sds)",
   max(abs(rp_s$scale - sd_res)) < 1e-8 && max(abs(rp_s$scale - sd_raw)) > 1e-2,
   sprintf("max|scale-residSD|=%.2e  max|scale-rawSD|=%.3f",
           max(abs(rp_s$scale - sd_res)), max(abs(rp_s$scale - sd_raw))))

## ---------------- T2: df guard uses DESIGN RANK (collinear age + age_copy)
meta12 <- meta[1:12, ]; meta12$sex <- factor(rep(c("F","M"), each = 6))
meta12$batch <- factor(rep(c("b1","b2"), each = 6))
x12 <- as_cell_program_data(list(NK = pb$NK[1:12, ]), sample_metadata = meta12)
g2a <- grab(fit_celltype_programs(x12, covariates = c("age","age_copy"), max_factors = 10,
                                  features = "variable", n_variable_genes = 100, backfit = FALSE))
g2b <- grab(fit_celltype_programs(x12, covariates = c("age","age_copy"), max_factors = 11,
                                  features = "variable", n_variable_genes = 100, backfit = FALSE))
ok("T2a collinear design: max_factors=10 does NOT trip guard (q_ct = rank-1 = 1, not ncol-1 = 2)",
   is.null(g2a$error) || !grepl("Too many factors", g2a$error),
   if (is.null(g2a$error)) "fit completed" else paste("error:", g2a$error))
ok("T2b max_factors=11 trips guard with q = 1 (rank-based)",
   !is.null(g2b$error) && grepl("Too many factors", g2b$error) && grepl("q = 1", g2b$error),
   if (!is.null(g2b$error)) g2b$error else "no error")

## ---------------- T3: residualization orthogonality (crossprod(R, X) ~ 0)
o1 <- max(abs(crossprod(R, Xres)))
Xres_c <- cd(ids, meta, c("age","age_copy"), "B")
R_c <- lm.fit(Xres_c, as.matrix(pb$B)[, fit_collin$preprocessing$feature_sets$B, drop = FALSE])$residuals
o2 <- max(abs(crossprod(R_c, Xres_c)))
ok("T3 residuals orthogonal to design (full-rank and collinear)",
   o1 < 1e-8 && o2 < 1e-8,
   sprintf("full-rank max|t(R)X| = %.2e ; collinear max|t(R)X| = %.2e", o1, o2))

## ---------------- T4: formula vs character equivalence
X_f <- cd(ids, meta, ~ age + sex + batch, "B")
X_c <- cd(ids, meta, c("age","sex","batch"), "B")
ok("T4 .cov_design formula == character", identical(X_f, X_c), sprintf("identical: %s", identical(X_f, X_c)))

## ---------------- T5: droplevels of unused factor levels
meta_lv <- meta; meta_lv$sex <- factor(as.character(meta_lv$sex), levels = c("F","M","X"))
X5 <- cd(ids, meta_lv, "sex", "B")
X5_ref <- stats::model.matrix(~ sex,
                              data = data.frame(sex = droplevels(meta_lv$sex[match(ids, meta_lv$observation_id)])))
ok("T5 unused level dropped, no all-zero dummy column (== independent sex-only design)",
   !any(colSums(abs(X5)) == 0) &&
     isTRUE(all.equal(unname(X5), unname(X5_ref), check.attributes = FALSE)) &&
     identical(colnames(X5), colnames(X5_ref)),
   sprintf("columns: %s", paste(colnames(X5), collapse = ", ")))

## ---------------- T6: constant-within-celltype covariates
g6  <- grab(X6 <- cd(ids, meta, c("site_const","age"), "B"))
g6b <- grab(X6b <- cd(ids, meta, "site_const", "B"))
ok("T6a constant covariate dropped with message",
   any(grepl("constant in cell type", g6$messages)) && !any(grepl("site_const", colnames(X6))),
   paste(g6$messages, collapse = " | "))
ok("T6b all-constant -> NULL design + 'not adjusted' message",
   is.null(X6b) && any(grepl("not adjusted", g6b$messages)), paste(g6b$messages, collapse = " | "))
meta_c8 <- data.frame(observation_id = ids[1:12], age = 50, sex = factor("F"), batch = factor("b1"))
x_c8 <- as_cell_program_data(list(CD8 = pb$CD8[1:12, ]), sample_metadata = meta_c8)
g6c <- grab(f8 <- fit_celltype_programs(x_c8, covariates = c("age","sex","batch"), max_factors = 3,
                                        features = "variable", n_variable_genes = 100, backfit = FALSE))
ok("T6c full fit with all-constant covariates proceeds unadjusted (coefficients NULL)",
   !is.null(f8) && is.null(f8$fits$CD8$preproc_replay$residualize_coefficients) &&
     any(grepl("not\\s+adjusted", g6c$messages)), paste(head(g6c$messages, 1), collapse = " | "))

## ---------------- T7: NA abort
meta_na <- meta; meta_na[5, "age"] <- NA_real_
g7 <- grab(cd(ids, meta_na, "age", "B"))
ok("T7 NA in covariate aborts", !is.null(g7$error) && grepl("Missing covariate values", g7$error),
   if (!is.null(g7$error)) g7$error else "no error")

## ---------------- T8: rank check nrow - rank >= 10
g8 <- grab(cd(ids[1:12], meta12, c("age","sex","batch"), "NK"))
ok("T8 df guard aborts at nrow - rank < 10 (12 donors, rank 4)",
   !is.null(g8$error) && grepl("Too few donors", g8$error),
   if (!is.null(g8$error)) g8$error else "no error")

## ---------------- T9: leverage warning for singleton level
meta_lvg <- meta; meta_lvg$batch <- factor(c(rep("b1",37), rep("b2",2), rep("b3",1)))
g9 <- grab(X9 <- cd(ids, meta_lvg, "batch", "B"))
hat <- rowSums(qr.Q(qr(X9))[, seq_len(qr(X9)$rank), drop = FALSE]^2)
ok("T9 leverage warning fires for singleton level (hat ~ 1)",
   any(grepl("leverage", g9$warnings)) && max(hat) > 0.99,
   sprintf("max hat = %.4f ; warning: %s", max(hat),
           if (length(g9$warnings)) head(g9$warnings,1) else "NONE"))

## ---------------- T10: fixed-mode design (intercept dropped, scaled, rank-pivoted)
Xc <- cd(ids, meta, c("age","sex","batch"), "B", intercept = FALSE)
ok("T10a fixed design: no intercept, centered, unit variance, full rank",
   !any(colnames(Xc) == "(Intercept)") && max(abs(colMeans(Xc))) < 1e-10 &&
     max(abs(apply(Xc, 2, sd) - 1)) < 1e-10 && qr(Xc)$rank == ncol(Xc),
   sprintf("columns: %s", paste(colnames(Xc), collapse = ", ")))
meta_s2 <- meta; meta_s2$sex2 <- meta_s2$sex
Xc2 <- cd(ids, meta_s2, c("sex","sex2"), "B", intercept = FALSE)
ok("T10b fixed design rank-pivots aliased columns (sex == sex2 -> 1 column)", ncol(Xc2) == 1L,
   sprintf("ncol = %d", ncol(Xc2)))

## ---------------- T11: .flash_fixed semantics
fl_B <- fit_fix$fits$B$flash
Xc_ref <- cd(ids, meta, ~ age + sex + batch, "B", intercept = FALSE)
ok("T11a fixed_L is exactly the fixed design",
   identical(unname(fl_B$fixed_L), unname(Xc_ref)) && identical(colnames(fl_B$fixed_L), colnames(Xc_ref)))
Yc <- scale(Yraw, center = TRUE, scale = FALSE)   # fit default is scale=FALSE; scale() defaults to scale=TRUE
d11 <- max(abs(Yc - fl_B$fixed_L %*% t(fl_B$fixed_F) - fit_fix$fits$B$Y_used))
ok("T11b Y_used == centered Y - fixed_L %*% t(fixed_F)  ('Y minus fitted covariate part')",
   d11 < 1e-8, sprintf("max|diff| = %.2e", d11))
F0_ols <- t(solve(crossprod(Xc_ref), crossprod(Xc_ref, Yc)))
cat(sprintf("       note: max|fixed_F - OLS init| = %.3f (fixed_F is the flash-estimated factor; OLS is only its initializer)\n",
            max(abs(fl_B$fixed_F - F0_ols))))
fit_none_k <- canonicalize_programs(fit_none)
fit_res_k  <- canonicalize_programs(fit_res)   # reused by T12 below
fit_fix_k  <- canonicalize_programs(fit_fix)   # reused by T18 below
K_free <- fl_B$n_factors
ok("T11c fixed columns NOT returned as free programs",
   ncol(fit_fix_k$scores$B) == K_free &&
     !any(grepl("age|sex|batch|Intercept", colnames(fit_fix_k$scores$B))),
   sprintf("K_free = %d, program names: %s", K_free,
           paste(head(colnames(fit_fix_k$scores$B)), collapse = ",")))
csn <- max(abs(apply(fit_none_k$scores$B, 2, function(z) suppressWarnings(cor(age, z)))))
csr <- max(abs(apply(fit_res_k$scores$B, 2, function(z) suppressWarnings(cor(age, z)))))
csf <- max(abs(apply(fit_fix_k$scores$B, 2, function(z) suppressWarnings(cor(age, z)))))
cat(sprintf("       age-vs-score |cor|: no-cov=%.3f  residualize=%.3f  fixed=%.3f\n", csn, csr, csf))
ok("T11d planted age module removed from free programs (both modes)", csr < 0.5*csn && csf < 0.5*csn)
selB <- fit_fix$preprocessing$feature_sets$B
planted <- selB %in% colnames(pb$B)[1:20]
if (sum(planted) > 5) {
  cfx <- cor(fl_B$fixed_F[planted, "age"], w_age[match(selB[planted], colnames(pb$B))])
  cat(sprintf("       cor(fixed_F age column, planted age loadings) = %.3f\n", cfx))
}
g11e <- grab(fit_celltype_programs(x0, covariates = ~ age + sex + batch, covariate_mode = "fixed",
                                   flash_control = list(tol = 1e-2), max_factors = 3,
                                   features = "variable", n_variable_genes = 100, backfit = FALSE))
ok("T11e fixed mode silently ignores flash_control (documented)", is.null(g11e$error),
   "doc: 'Requires center = TRUE, ignores flash_control'")

## ---------------- T12: projection replay reproduces in-sample adjustment
cellsB  <- setNames(rep("B", 40), ids)
donorsB <- setNames(ids, ids)
g12 <- grab(pr <- project_programs_to_cells(fit_res_k, pb$B, cell_types = cellsB,
                                            donors = donorsB, sample_metadata = meta))
ok("T12a projection of training data runs", is.null(g12$error),
   if (!is.null(g12$error)) g12$error else "")
shared <- intersect(fit_res_k$fits$B$genes, colnames(pb$B))
Xr <- stats::model.matrix(~ age + sex + batch, meta[match(ids, meta$observation_id), c("age","sex","batch")])
beta <- fit_res_k$fits$B$preproc_replay$residualize_coefficients
Q <- as.matrix(pb$B)[ids, shared, drop = FALSE] -
  Xr %*% beta[fit_res_k$fits$B$preproc_replay$design_columns, shared, drop = FALSE]
idx <- match(shared, fit_res_k$fits$B$genes)
Q <- sweep(Q, 2, fit_res_k$fits$B$preproc_replay$center[idx], "-")
Q <- sweep(Q, 2, fit_res_k$fits$B$preproc_replay$scale[idx], "/")
d12 <- max(abs(Q - fit_res_k$fits$B$Y_used[, shared]))
ok("T12b replayed adjustment == in-sample Y_used (same donors)", d12 < 1e-8,
   sprintf("max|diff| = %.2e", d12))
ScB <- fit_res_k$scores$B
pcols <- paste0("B_", seq_len(ncol(ScB)))   # activity spans ALL cell types' programs; compare B's
if (all(pcols %in% colnames(pr$activity))) {
  cors <- sapply(seq_len(ncol(ScB)), function(k) cor(pr$activity[, pcols[k]], ScB[, k]))
  ok("T12c projected activities match donor scores (B program columns)",
     min(abs(cors)) > 0.9,
     sprintf("min|cor| = %.4f (activity matrix has %d cols across all cell types)",
             min(abs(cors)), ncol(pr$activity)))
} else {
  ok("T12c projected activities match donor scores (B program columns)", FALSE,
     "B program columns missing from activity matrix")
}
ord <- sample(40)
prS <- project_programs_to_cells(fit_res_k, pb$B[ord, ], cell_types = cellsB[ord],
                                 donors = donorsB[ord], sample_metadata = meta)
d12d <- max(abs(prS$activity[ids, pcols] - pr$activity[ids, pcols]))
ok("T12d replay order-independent", d12d < 1e-8, sprintf("max|diff| = %.2e", d12d))

## ---------------- T13: design-mismatch fail-closed errors
meta_bad <- meta; meta_bad$batch <- factor(ifelse(seq_len(40) > 30, "b9", as.character(meta$batch)))
g13a <- grab(project_programs_to_cells(fit_res_k, pb$B, cell_types = cellsB,
                                       donors = donorsB, sample_metadata = meta_bad))
ok("T13a new covariate level among query donors aborts",
   !is.null(g13a$error) && grepl("does not match the design", g13a$error),
   if (!is.null(g13a$error)) g13a$error else "NO ERROR (silent!)")
meta_lost <- meta; meta_lost$batch <- factor(rep("b1", 40))
g13b <- grab(project_programs_to_cells(fit_res_k, pb$B, cell_types = cellsB,
                                       donors = donorsB, sample_metadata = meta_lost))
ok("T13b lost covariate level among query donors aborts (fail-closed; obscure msg -> P3 note)",
   !is.null(g13b$error),
   sprintf("error: %s", substr(if (!is.null(g13b$error)) g13b$error else "NO ERROR (silent!)", 1, 90)))
g13c <- grab(project_programs_to_cells(fit_res_k, pb$B, cell_types = cellsB,
                                       donors = donorsB, sample_metadata = meta[1:39, ]))
ok("T13c query donor missing from sample_metadata aborts",
   !is.null(g13c$error) && grepl("missing from", g13c$error),
   if (!is.null(g13c$error)) g13c$error else "NO ERROR")
g13d <- grab(project_programs_to_cells(fit_res_k, pb$B, cell_types = cellsB))
ok("T13d missing donors/sample_metadata aborts when covariates were residualized",
   !is.null(g13d$error) && grepl("required to replay", g13d$error),
   if (!is.null(g13d$error)) g13d$error else "NO ERROR")
meta_nocol <- meta[, c("observation_id","age","sex")]
g13e <- grab(project_programs_to_cells(fit_res_k, pb$B, cell_types = cellsB,
                                       donors = donorsB, sample_metadata = meta_nocol))
ok("T13e metadata missing a covariate column fails (fail-closed)",
   !is.null(g13e$error),
   sprintf("error class: %s", if (!is.null(g13e$error)) substr(g13e$error, 1, 60) else "NONE"))

## ---------------- T14: preprocessing_spec interplay
spec <- freeze_preprocessing(x0, features = "variable", n_variable_genes = 100,
                             covariates = ~ age + sex + batch)
ok("T14a formula stored deparsed", is.character(spec$covariates) && grepl("^\\s*~", spec$covariates),
   paste(spec$covariates))
g14b <- grab(fspec <- fit_celltype_programs(x0, max_factors = 3, backfit = FALSE,
                                            preprocessing_spec = spec))
ok("T14b spec refit replays covariates (rehydrated formula) + records digest",
   is.null(g14b$error) && inherits(fspec$preprocessing$covariates, "formula") &&
     identical(fspec$preprocessing$frozen_spec, spec$digest),
   if (!is.null(g14b$error)) g14b$error else "")
g14c <- grab(fit_celltype_programs(x0, covariates = ~ age, max_factors = 3, preprocessing_spec = spec))
ok("T14c explicit covariates + frozen spec -> error", !is.null(g14c$error) && grepl("frozen", g14c$error),
   if (!is.null(g14c$error)) g14c$error else "NO ERROR (silent drift!)")
g14d <- grab(fit_celltype_programs(x0, covariate_mode = "fixed", max_factors = 3, preprocessing_spec = spec))
ok("T14d explicit covariate_mode + frozen spec -> error",
   !is.null(g14d$error) && grepl("frozen", g14d$error),
   if (!is.null(g14d$error)) g14d$error else "NO ERROR (silent drift!)")
Xspec <- cd(ids, meta, cellprograms:::.spec_covariates(spec$covariates), "B")
ok("T14e spec-rehydrated formula design == direct formula design", identical(Xspec, X_f))
spec_bad <- freeze_preprocessing(x0, features = "variable", n_variable_genes = 100,
                                 covariates = ~ age + sex, covariate_mode = "fixed", center = FALSE)
g14f <- grab(fit_celltype_programs(x0, max_factors = 3, backfit = FALSE, preprocessing_spec = spec_bad))
g14f_defect <- is.null(g14f$error)
ok("T14f [DEFECT MARKER P2-3] frozen spec bypasses fixed-mode center=TRUE guard (passes while defect present; invert when fixed)",
   g14f_defect,
   if (g14f_defect) "confirmed: fixed mode ran with center=FALSE via frozen spec" else "guard now fires (defect fixed)")

## ---------------- T15: per-cell-type design rebuild
pb2 <- list(B = pb$B, CD4 = pb$CD4[1:30, ])
x2 <- as_cell_program_data(pb2, sample_metadata = meta)
g15 <- grab(f2 <- fit_celltype_programs(x2, covariates = c("site","sex"), max_factors = 3,
                                        features = "variable", n_variable_genes = 100, backfit = FALSE))
dcB <- f2$fits$B$preproc_replay$design_columns
dcC <- f2$fits$CD4$preproc_replay$design_columns
ok("T15 site constant in CD4 dropped there, kept in B (own designs)",
   any(grepl("constant in cell type \u00abCD4\u00bb|constant in cell type .CD4", g15$messages)) &&
     !any(grepl("site", dcC)) && any(grepl("siteS2", dcB)),
   sprintf("CD4 columns: %s | B columns: %s", paste(dcC, collapse = ","), paste(dcB, collapse = ",")))

## ---------------- T16: collinear design -> NA coefficients -> projection replay
beta16 <- fit_collin$fits$B$preproc_replay$residualize_coefficients
ok("T16a collinear design stores NA coefficient for aliased column",
   anyNA(beta16["age_copy", ]), sprintf("NA rows: %s",
   paste(rownames(beta16)[apply(is.na(beta16), 1, any)], collapse = ",")))
fit_collin_k <- canonicalize_programs(fit_collin)
g16 <- grab(pr16 <- project_programs_to_cells(fit_collin_k, pb$B, cell_types = cellsB,
                                              donors = donorsB, sample_metadata = meta))
if (!is.null(g16$error)) {
  cat(sprintf("       T16b outcome: projection ERRORS: %s\n", substr(g16$error, 1, 120)))
  ok("T16b collinear-design projection fails closed (error, not silent NA)", TRUE, "error raised")
} else {
  n_na <- sum(is.na(pr16$activity))
  cat(sprintf("       T16b outcome: projection returned %d NA activities out of %d\n",
              n_na, length(pr16$activity)))
  ok("T16b collinear-design projection silently degraded", FALSE,
     "NO error raised -- inspect NA count above")
}

## ---------------- T17: replay center/scale misalignment when zero-variance genes are dropped
pb3 <- pb; pb3$B[, "G1"] <- 5; pb3$B[, "G2"] <- pb3$B[, "G2"] + 10
x3 <- as_cell_program_data(pb3, sample_metadata = meta)
g17 <- grab(f3 <- fit_celltype_programs(x3, covariates = ~ age + sex, max_factors = 3,
                                        features = "all", backfit = FALSE, feature_scale = "logcounts"))
rp3 <- f3$fits$B$preproc_replay
ok("T17a replay$center length (pre-drop) != nrow genes (post-drop)",
   length(rp3$center) == 150L && length(f3$fits$B$genes) == 149L,
   sprintf("length(center) = %d, ngenes = %d, drop message: %s", length(rp3$center),
           length(f3$fits$B$genes), any(grepl("zero-variance", g17$messages))))
# no-covariate fit, same data: default path
g17b <- grab(f4 <- fit_celltype_programs(x3, max_factors = 3, features = "all",
                                         backfit = FALSE, feature_scale = "logcounts"))
rp4 <- f4$fits$B$preproc_replay
sh4 <- intersect(f4$fits$B$genes, colnames(pb3$B))
Qc4 <- as.matrix(pb3$B)[ids, sh4, drop = FALSE]
i4 <- match(sh4, f4$fits$B$genes)
Qc4 <- sweep(Qc4, 2, rp4$center[i4], "-"); Qc4 <- sweep(Qc4, 2, rp4$scale[i4], "/")
Qo4 <- as.matrix(pb3$B)[ids, sh4, drop = FALSE]
Qo4 <- sweep(Qo4, 2, rp4$center[-1], "-"); Qo4 <- sweep(Qo4, 2, rp4$scale[-1], "/")
dc4 <- max(abs(Qc4 - f4$fits$B$Y_used[, sh4])); do4 <- max(abs(Qo4 - f4$fits$B$Y_used[, sh4]))
ok("T17b DEFECT: code-replayed centering misaligned vs Y_used; aligned replay matches",
   dc4 > 1 && do4 < 1e-8,
   sprintf("code-replay max|diff| = %.3f ; correctly-aligned replay max|diff| = %.2e", dc4, do4))
# covariate + scale=TRUE + constant gene: does it crash?
g17c <- grab(fit_celltype_programs(x3, covariates = ~ age + sex, scale = TRUE, max_factors = 3,
                                   features = "all", backfit = FALSE))
cat(sprintf("       T17c covariates+scale=TRUE+constant gene: %s\n",
            if (!is.null(g17c$error)) paste("ERROR:", substr(g17c$error, 1, 100)) else "completed"))

## ---------------- T18: fixed-covariate fits refuse projection (documented)
g18 <- grab(project_programs_to_cells(fit_fix_k, pb$B, cell_types = cellsB))
ok("T18 fixed-covariate fit refuses projection",
   !is.null(g18$error) && grepl("fixed-covariate|manifest", g18$error),
   if (!is.null(g18$error)) substr(g18$error, 1, 100) else "NO ERROR")

## ---------------- T19: stability refits carry covariates
g19 <- grab(st <- assess_program_stability(fit_res_k, n_boot = 2, resample_frac = 0.8, progress = FALSE))
ok("T19 stability subsample refits run with covariates",
   is.null(g19$error) && !is.null(st$stability) && nrow(st$stability$per_program) > 0,
   if (!is.null(g19$error)) g19$error else sprintf("%d programs tracked", nrow(st$stability$per_program)))

## ---------------- T20: union_panel_loadings replay
g20a <- grab(up_r <- union_panel_loadings(fit_res_k))
g20b <- grab(up_f <- union_panel_loadings(fit_fix_k))
fin <- function(m) { if (is.list(m)) all(vapply(m, function(mm) all(is.finite(mm) | is.na(mm)), TRUE)) else all(is.finite(m) | is.na(m)) }
ok("T20 union_panel_loadings replays residualize and fixed adjustment",
   is.null(g20a$error) && is.null(g20b$error) && fin(up_r$loadings) && fin(up_f$loadings),
   sprintf("errors: %s / %s", if (!is.null(g20a$error)) g20a$error else "none",
           if (!is.null(g20b$error)) g20b$error else "none"))

## ---------------- T21: formula interaction involving a within-ct-constant covariate
meta12c <- meta12; meta12c$sex <- factor(rep("F", 12))
g21 <- grab(cd(ids[1:12], meta12c, ~ age + sex + age:sex, "NK"))
ok("T21 [DEFECT MARKER P3] interaction with within-ct-constant covariate errors obscurely (passes while defect present; invert when fixed)",
   !is.null(g21$error),
   if (!is.null(g21$error)) paste("confirmed ERROR:", substr(g21$error, 1, 100)) else "handled gracefully (defect fixed?)")

## ---------------- T22: character vs formula full-fit equivalence
ok("T22 character and formula covariates give identical Y_used and design columns",
   identical(fit_res$fits$B$Y_used, fit_res_c$fits$B$Y_used) &&
     identical(fit_res$fits$B$preproc_replay$design_columns,
               fit_res_c$fits$B$preproc_replay$design_columns))

cat(sprintf("\n== SUMMARY: %d PASS, %d FAIL ==\n", PASS, FAIL))
cat("\n-- sessionInfo --\n")
print(utils::sessionInfo()$running)
cat("cellprograms", as.character(packageVersion("cellprograms")), "\n")
