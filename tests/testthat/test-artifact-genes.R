# Toy fit with a planted V(D)J/genotype-marker program; loadings are planted
# directly, no flashier needed. A_1 carries IGHV1-2/IGKC/GSTM1, B_2 carries
# XIST/RPS4Y1, everything else is small random noise.
.artifact_toy <- function(seed = 5) {
  set.seed(seed)
  genes <- c(paste0("g", 1:36), "IGHV1-2", "IGKC", "GSTM1", "XIST", "RPS4Y1")
  artifacts <- c("IGHV1-2", "IGKC", "GSTM1", "XIST", "RPS4Y1")
  mk <- function(ct, plant) {
    W <- matrix(rnorm(length(genes) * 2, sd = 0.1), length(genes), 2,
                dimnames = list(genes, paste0(ct, "_", 1:2)))
    W[artifacts, ] <- 0.001                # artifacts never reach the top unless planted
    if (identical(plant, "vdj")) {
      W[c("IGHV1-2", "IGKC", "GSTM1"), 1] <- c(5, 4, 3)
    } else if (identical(plant, "sex")) {
      W[c("XIST", "RPS4Y1"), 2] <- c(5, 4)
    }
    W
  }
  structure(list(
    loadings = list(A = mk("A", "vdj"), B = mk("B", "sex")),
    cell_types = c("A", "B")
  ), class = c("cell_program_fit", "cell_program_data", "list"))
}

test_that("artifact_gene_sets stores the documented sets and .expand_vdj expands the regex", {
  expect_named(artifact_gene_sets, c("vdj", "genotype_marker", "sex"))
  expect_length(artifact_gene_sets$vdj, 1L)          # pattern, not a gene list
  expect_setequal(artifact_gene_sets$genotype_marker,
                  c("GSTM1", "GSTT1", "GSTT2", "UTS2", "ERAP2", "CCL3L1", "CCL3L3",
                    "CCL4L1", "CCL4L2", "HLA-DRB3", "HLA-DRB4", "HLA-DRB5",
                    "HLA-DQA1", "HLA-DQA2", "HLA-DQB1", "HLA-DQB2"))
  expect_setequal(artifact_gene_sets$sex,
                  c("XIST", "RPS4Y1", "DDX3Y", "UTY", "KDM5D", "EIF1AY",
                    "NLGN4Y", "TMSB4Y", "USP9Y", "ZFY"))
  universe <- c("IGHV1-2", "IGKC", "IGHG1", "JCHAIN", "TRBV7-2", "GSTM1", "MALAT1")
  expect_setequal(.expand_vdj(universe), c("IGHV1-2", "IGKC", "IGHG1", "JCHAIN", "TRBV7-2"))
  expect_setequal(artifact_genes(universe, "vdj"), .expand_vdj(universe))
  expect_setequal(artifact_genes(universe), c("IGHV1-2", "IGKC", "IGHG1", "JCHAIN", "TRBV7-2", "GSTM1"))
  expect_error(artifact_genes(universe, "nope"), "Unknown artifact gene set")
})

test_that("flag_artifact_genes flags a planted V(D)J-carrying program (cmp = NULL)", {
  fit <- .artifact_toy()
  fl <- flag_artifact_genes(fit)
  expect_equal(nrow(fl), 4L)                          # one row per program
  expect_setequal(names(fl), c("cell_type", "program_id",
                               "vdj_frac", "vdj_flag", "genotype_marker_frac", "genotype_marker_flag"))
  a1 <- fl[fl$program_id == "A_1", ]
  expect_true(a1$vdj_flag)
  expect_true(a1$genotype_marker_flag)
  expect_equal(a1$vdj_frac, 2 / 12)                   # IGHV1-2, IGKC in top 12
  expect_equal(a1$genotype_marker_frac, 1 / 12)       # GSTM1
  b1 <- fl[fl$program_id == "B_1", ]
  expect_false(b1$vdj_flag)
  expect_false(b1$genotype_marker_flag)
  expect_equal(b1$vdj_frac, 0)
  # sex set is opt-in via `sets`
  fls <- flag_artifact_genes(fit, sets = "sex")
  expect_true(fls$sex_flag[fls$program_id == "B_2"])
  expect_false(fls$sex_flag[fls$program_id == "A_1"])
  expect_error(flag_artifact_genes(fit, sets = "nope"), "Unknown artifact gene set")
})

test_that("flag_artifact_genes pools member top genes per CMP", {
  fit <- .artifact_toy()
  tab <- data.frame(cmp = c(1L, 1L, 2L),
                    cell_type = c("A", "B", "B"),
                    program_id = c("A_1", "B_1", "B_2"),
                    stringsAsFactors = FALSE)
  fl <- flag_artifact_genes(fit, tab, sets = c("vdj", "genotype_marker", "sex"))
  expect_equal(nrow(fl), 2L)                          # one row per CMP
  expect_equal(fl$n_members, c(2L, 1L))
  expect_true(fl$vdj_flag[fl$cmp == 1])               # carried by member A_1
  expect_true(fl$genotype_marker_flag[fl$cmp == 1])
  expect_false(fl$sex_flag[fl$cmp == 1])
  expect_true(fl$sex_flag[fl$cmp == 2])
  expect_false(fl$vdj_flag[fl$cmp == 2])
  # empty CMP table -> zero-row result with the flag columns
  e <- flag_artifact_genes(fit, tab[0, ], sets = "vdj")
  expect_equal(nrow(e), 0L)
  expect_true(all(c("vdj_frac", "vdj_flag") %in% names(e)))
})

test_that("exclude_genes is applied before variable-gene selection", {
  set.seed(2)
  n <- 20; g <- 30
  Y <- matrix(rnorm(n * g), n, g, dimnames = list(paste0("s", 1:n), paste0("g", 1:g)))
  Y[, "g1"] <- Y[, "g1"] * 10                         # g1 is the most variable gene
  cpd <- as_cell_program_data(list(CT1 = Y))
  fs0 <- .select_features(cpd, "variable", 10)
  expect_true("g1" %in% fs0$CT1)                      # baseline: g1 takes an HVG slot
  fs1 <- .select_features(cpd, "variable", 10, exclude_genes = "g1")
  expect_false("g1" %in% fs1$CT1)
  expect_length(fs1$CT1, 10L)                         # slot is reallocated, not lost
  v <- apply(Y[, setdiff(colnames(Y), "g1"), drop = FALSE], 2, stats::var)
  expect_setequal(fs1$CT1, names(sort(v, decreasing = TRUE))[1:10])
  # "all" and explicit-list modes drop excluded genes too
  expect_false("g1" %in% .select_features(cpd, "all", exclude_genes = "g1")$CT1)
  expl <- .select_features(cpd, list(CT1 = colnames(Y)[1:15]), exclude_genes = "g1")
  expect_setequal(expl$CT1, colnames(Y)[2:15])
  # exclude_genes = NULL is a no-op
  expect_identical(.select_features(cpd, "variable", 10),
                   .select_features(cpd, "variable", 10, exclude_genes = NULL))
})

test_that("exclude_genes removes genes from the fitted panel; default unchanged", {
  skip_if_not_installed("flashier")
  skip_if_not_installed("ebnm")
  set.seed(4)
  n <- 24; g <- 40
  Y <- matrix(rnorm(n * g, sd = 0.5), n, g, dimnames = list(paste0("s", 1:n), paste0("g", 1:g)))
  colnames(Y)[1:3] <- c("IGHV1-2", "IGKC", "GSTM1")
  Y[, 1:3] <- Y[, 1:3] * 6 + 3                        # artifact genes are highly variable
  z <- rnorm(n)
  Y[, 10:19] <- Y[, 10:19] + outer(z, rep(1, 10)) * 2 # planted program
  cpd <- as_cell_program_data(list(CT1 = Y))
  excl <- artifact_genes(colnames(Y), c("vdj", "genotype_marker"))
  expect_setequal(excl, c("IGHV1-2", "IGKC", "GSTM1"))
  fit <- canonicalize_programs(fit_celltype_programs(
    cpd, features = "variable", n_variable_genes = 25, exclude_genes = excl, seed = 1))
  expect_false(any(excl %in% fit$fits$CT1$genes))
  expect_false(any(excl %in% rownames(fit$loadings$CT1)))
  expect_equal(length(fit$fits$CT1$genes), 25L)       # excluded genes freed their HVG slots
  expect_equal(fit$preprocessing$exclude_genes, excl)
  # default (exclude_genes = NULL): artifact genes stay in the panel
  fit0 <- canonicalize_programs(fit_celltype_programs(
    cpd, features = "variable", n_variable_genes = 25, seed = 1))
  expect_true(all(excl %in% fit0$fits$CT1$genes))
  expect_null(fit0$preprocessing$exclude_genes)
  expect_error(fit_celltype_programs(cpd, exclude_genes = 1:3), "exclude_genes")
})
