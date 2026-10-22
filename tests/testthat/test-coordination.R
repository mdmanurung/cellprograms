# Toy canonicalized fit: no flashier needed; scores/loadings planted directly.
.cmp_toy <- function(n = 80, k = 3, g = 200, seed = 1, couple = TRUE, site_shift = 0, context_only = FALSE,
                     cts = c("A", "B", "C")) {
  set.seed(seed)
  ids <- sprintf("d%03d", seq_len(n))
  site <- rep(c("s1", "s2"), length.out = n)
  ctx <- rep(c(0, 1), each = n / 2)
  z <- rnorm(n)
  genes <- sprintf("G%03d", seq_len(g))
  scores <- loadings <- pb <- list()
  for (ct in cts) {
    S <- matrix(rnorm(n * k), n, k, dimnames = list(ids, paste0(ct, "_", seq_len(k))))
    S <- S + site_shift * (site == "s1")
    if (couple && ct %in% cts[1:3]) {
      j <- match(ct, cts)                       # coupled program index differs by cell type
      S[, j] <- (if (context_only) 2 * ctx else z) * (if (ct == cts[3]) -1 else 1) + rnorm(n, sd = 0.6)
    }
    W <- matrix(rnorm(g * k, sd = 0.1), g, k, dimnames = list(genes, colnames(S)))
    W[(1:20) + 20 * (match(ct, cts) - 1), ] <- W[(1:20) + 20 * (match(ct, cts) - 1), ] + 1   # distinct genes
    scores[[ct]] <- S; loadings[[ct]] <- W
    pb[[ct]] <- matrix(rnorm(n * g, 5), n, g, dimnames = list(ids, genes))
  }
  structure(list(scores = scores, loadings = loadings, cell_types = cts, observation_ids = ids, pseudobulk = pb,
                 sample_metadata = data.frame(observation_id = ids, site = site, ctx = ctx, stringsAsFactors = FALSE)),
            class = c("cell_program_fit", "cell_program_data", "list"))
}

test_that("planted cross-cell-type coupling with distinct genes is recovered as one CMP", {
  fit <- .cmp_toy()
  res <- find_coordinated_programs(fit, strata = "site", n_perm = 99, ambient = "none")
  expect_s3_class(res, "cmp_result")
  expect_equal(length(unique(res$cmps$cmp)), 1L)
  expect_setequal(paste(res$cmps$cell_type, res$cmps$program_id), c("A A_1", "B B_2", "C C_3"))
  expect_equal(res$cmps$sign[res$cmps$cell_type == "C"], -1)
  expect_true(all(res$cmps$node_p < 0.05))
  expect_equal(cmp_gene_overlap(fit, res)$class, "distinct-genes")
  sc <- cmp_scores(fit, res)
  expect_equal(dim(sc), c(80L, 1L))
  expect_gt(abs(stats::cor(sc[, 1], fit$scores$A[, 1])), 0.7)
})

test_that("no CMPs are called on uncoupled data (FWER)", {
  calls <- vapply(1:15, function(s) {
    fit <- .cmp_toy(seed = 100 + s, couple = FALSE)
    nrow(find_coordinated_programs(fit, strata = "site", n_perm = 99, ambient = "none", seed = s)$cmps) > 0
  }, TRUE)
  expect_lte(sum(calls), 3L)
})

test_that("strata remove site-driven pseudo-coordination", {
  fit <- .cmp_toy(couple = FALSE, site_shift = 2, seed = 7)
  with_strata <- find_coordinated_programs(fit, strata = "site", n_perm = 99, ambient = "none")
  without <- find_coordinated_programs(fit, n_perm = 99, ambient = "none")
  expect_equal(nrow(with_strata$cmps), 0L)
  expect_gt(nrow(without$cmps), 0L)
})

test_that("within mode removes coordination explained by context", {
  fit <- .cmp_toy(context_only = TRUE, seed = 3)
  marg <- find_coordinated_programs(fit, strata = "site", context = "ctx", n_perm = 99, ambient = "none")
  within <- find_coordinated_programs(fit, strata = "site", context = "ctx", mode = "within", n_perm = 99,
                                      ambient = "none")
  expect_gt(nrow(marg$cmps), 0L)
  expect_equal(nrow(within$cmps), 0L)
  expect_error(find_coordinated_programs(fit, mode = "within"), "requires `context`")
})

test_that("ambient-like programs are dropped or flagged", {
  fit <- .cmp_toy(seed = 11)
  # make A_1 (a coupled program) track B's expression specificity
  spec <- rnorm(200); names(spec) <- rownames(fit$loadings$A)
  prof <- list(A = setNames(rep(5, 200), names(spec)), B = 5 + spec, C = setNames(rep(5, 200), names(spec)))
  fit$loadings$A[, "A_1"] <- spec + rnorm(200, sd = 0.2)
  amb <- ambient_program_scores(fit, prof)
  expect_gt(amb$ambient_score[amb$program_id == "A_1"], 0.8)
  expect_equal(amb$source[amb$program_id == "A_1"], "B")
  dropped <- find_coordinated_programs(fit, strata = "site", n_perm = 99, profiles = prof)
  expect_false("A_1" %in% dropped$cmps$program_id)
  flagged <- find_coordinated_programs(fit, strata = "site", n_perm = 99, profiles = prof, ambient = "flag")
  expect_true(flagged$cmps$ambient_flag[flagged$cmps$program_id == "A_1"])
})

test_that("LR, theme and context annotations behave", {
  fit <- .cmp_toy(seed = 5, g = 600)
  res <- find_coordinated_programs(fit, strata = "site", n_perm = 99, ambient = "none")
  # plant a ligand at the top of A_1 and its receptor at the top of B_2
  fit$loadings$A["G150", "A_1"] <- 5; fit$loadings$B["G151", "B_2"] <- 5
  # sparse resource (as in real LR databases): the label-permutation null rarely re-creates the planted pair
  set.seed(2)
  lr <- data.frame(ligand = c("G150", sample(sprintf("G%03d", setdiff(61:349, 150:151)), 300, TRUE)),
                   receptor = c("G151", sample(sprintf("G%03d", 352:600), 300, TRUE)))
  hits <- annotate_cmp_lr(fit, res, lr, n_perm = 199)
  expect_true(any(hits$ligand == "G150" & hits$receptor == "G151"))
  expect_true(all(c("donor_cor", "n_donors") %in% names(hits)))
  expect_equal(nrow(annotate_cmp_lr(fit, res, lr[-1, ], n_perm = 199)), 0L)
  th <- annotate_cmp_themes(fit, res, list(setA = sprintf("G%03d", 1:20), setNone = sprintf("G%03d", 181:200)))
  expect_equal(th$set[th$cell_type == "A"], "setA")
  ca <- cmp_context_association(fit, res, "ctx")
  expect_equal(nrow(ca), 1L)
  expect_true(all(c("estimate", "p") %in% names(ca)))
})

test_that("input validation", {
  fit <- .cmp_toy()
  expect_error(find_coordinated_programs(fit, strata = "nope"), "not found")
  expect_error(find_coordinated_programs(fit, n_perm = 5), "too small")
  bad <- fit; bad$scores <- NULL
  expect_error(find_coordinated_programs(bad), "canonicalize_programs")
  expect_output(print(find_coordinated_programs(fit, strata = "site", n_perm = 99, ambient = "none")), "CMP 1")
})

test_that("LR annotation with FDR control recovers the planted pair and returns BH-adjusted p", {
  fit <- .cmp_toy(seed = 5, g = 600)
  res <- find_coordinated_programs(fit, strata = "site", n_perm = 99, ambient = "none")
  fit$loadings$A["G150", "A_1"] <- 5; fit$loadings$B["G151", "B_2"] <- 5
  set.seed(2)
  lr <- data.frame(ligand = c("G150", sample(sprintf("G%03d", setdiff(61:349, 150:151)), 300, TRUE)),
                   receptor = c("G151", sample(sprintf("G%03d", 352:600), 300, TRUE)))
  hf <- annotate_cmp_lr(fit, res, lr, n_perm = 199, control = "fdr")
  expect_true(any(hf$ligand == "G150" & hf$receptor == "G151"))
  expect_true(all(hf$p >= 0 & hf$p < 0.05))
  expect_equal(nrow(annotate_cmp_lr(fit, res, lr[-1, ], n_perm = 199, control = "fdr")), 0L)
  expect_error(annotate_cmp_lr(fit, res, lr, control = "bonferroni"))
})

test_that("context association drops donors with missing context and constant covariates", {
  fit <- .cmp_toy(seed = 5)
  res <- find_coordinated_programs(fit, strata = "site", n_perm = 99, ambient = "none")
  full <- cmp_context_association(fit, res, "ctx")
  fit$sample_metadata$ctx[1:10] <- NA
  ca <- cmp_context_association(fit, res, "ctx")
  expect_equal(ca$n, full$n - 10L)
  expect_true(is.finite(ca$p))
  fit$sample_metadata$const <- "x"
  cc <- cmp_context_association(fit, res, "ctx", covariates = c("site", "const"))
  expect_equal(nrow(cc), nrow(ca))
  expect_true(is.finite(cc$p))
})

test_that("a one-time message suggests max_factors = 30 only when a cell type hit a lower cap", {
  old <- options(cellprograms.cmp_cap_message_shown = NULL); on.exit(options(old))
  fit <- .cmp_toy(seed = 1)
  fit$preprocessing <- list(max_factors = 5)                  # toy has 3 programs per cell type: cap not reached
  expect_no_message(find_coordinated_programs(fit, strata = "site", n_perm = 99, ambient = "none"))
  fit$preprocessing <- list(max_factors = 3)                  # cap reached
  expect_message(find_coordinated_programs(fit, strata = "site", n_perm = 99, ambient = "none"), "max_factors = 30")
  expect_no_message(find_coordinated_programs(fit, strata = "site", n_perm = 99, ambient = "none"))   # once per session
  options(cellprograms.cmp_cap_message_shown = NULL)
  fit$preprocessing <- list(max_factors = 30)
  expect_no_message(find_coordinated_programs(fit, strata = "site", n_perm = 99, ambient = "none"))
})

test_that("empty score set after min_donors/ambient filtering fails closed (was: 'argument is of length zero')", {
  fit <- .cmp_toy(n = 6)   # 6 observed donors < min_donors = 20 default
  expect_error(
    find_coordinated_programs(fit, n_perm = 99, ambient = "none"),
    "No cell type retained programs"
  )
  expect_error(
    find_coordinated_programs(fit, n_perm = 99, ambient = "none"),
    "min_donors = 20"
  )
  # the same guard covers the ambient-drop path emptying every cell type
  expect_error(
    find_coordinated_programs(fit, n_perm = 99, ambient = "drop", ambient_threshold = -1),
    "No cell type retained programs"
  )
})
