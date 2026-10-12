# Uses .cmp_toy() from test-coordination.R (testthat sources test files in one environment per file,
# so the helper is redefined here via sys.source of that file's first definition if absent).
if (!exists(".cmp_toy")) {
  ex <- parse(test_path("test-coordination.R"))
  eval(ex[[which(vapply(ex, function(e) identical(as.character(e[[2]])[1], ".cmp_toy"), TRUE))[1]]])
}

.toy_net <- function() {
  g <- function(i) sprintf("G%03d", i)
  rbind(data.frame(source = "TFA", target = g(1:20), mor = 1),
        data.frame(source = "TFS", target = g(c(1:10, 21:30)), mor = 1),
        data.frame(source = "TFC", target = g(41:60), mor = 1),
        data.frame(source = "TFN", target = g(120:160), mor = 1),
        data.frame(source = "TFI", target = g(1:20), mor = -1))
}

test_that("prepare_collectri merges complexes and signs inhibition", {
  raw <- data.frame(source = c("P1", "COMPLEX:P2_P3", "COMPLEX:P4_P5", "P6", "P6"),
                    source_genesymbol = c("STAT1", "JUN_FOS", "RELA_NFKB1", "MYC", "MYC"),
                    target_genesymbol = c("IRF1", "MMP1", "NFKBIA", "CDKN1A", "CDKN1A"),
                    consensus_stimulation = c(TRUE, TRUE, TRUE, FALSE, FALSE),
                    consensus_inhibition = c(FALSE, FALSE, FALSE, TRUE, TRUE))
  net <- prepare_collectri(raw)
  expect_equal(net$source, c("STAT1", "AP1", "NFKB", "MYC"))
  expect_equal(net$mor, c(1, 1, 1, -1))
})

test_that("annotate_cmp_tf recovers planted regulons with CMP-oriented signs", {
  fit <- .cmp_toy()
  res <- find_coordinated_programs(fit, strata = "site", n_perm = 99, ambient = "none")
  tf <- annotate_cmp_tf(fit, res, .toy_net(), n_perm = 200)
  expect_s3_class(tf, "cmp_tf")
  m <- tf$members
  get <- function(ct, f) m[m$cell_type == ct & m$tf == f, ]
  expect_true(get("A", "TFA")$t > 0 && get("A", "TFA")$sig)
  expect_true(get("A", "TFI")$t < 0 && get("A", "TFI")$sig)
  expect_true(get("C", "TFC")$t < 0 && get("C", "TFC")$sig)   # C enters the CMP with sign -1
  expect_false(any(m$sig[m$tf == "TFN"]))   # depleted vs a skewed background, but not coherently signed
  s <- tf$summary
  expect_equal(s$class[s$tf == "TFS"], "shared")
  expect_equal(s$class[s$tf == "TFA"], "specific")
  expect_false("TFN" %in% s$tf)
})

test_that("annotate_cmp_tf p-values are roughly uniform for unrelated regulons", {
  fit <- .cmp_toy(seed = 11)
  res <- find_coordinated_programs(fit, strata = "site", n_perm = 99, ambient = "none")
  set.seed(5)
  net <- do.call(rbind, lapply(1:40, function(j)
    data.frame(source = paste0("R", j), target = sample(sprintf("G%03d", 61:200), 15), mor = 1)))
  tf <- annotate_cmp_tf(fit, res, net, n_perm = 200)
  expect_lt(mean(tf$members$p < 0.05), 0.15)
  expect_lte(sum(tf$members$sig), 2L)
})

test_that("cmp_triage flags nuisance genes, replication and context", {
  fit <- .cmp_toy()
  res <- find_coordinated_programs(fit, strata = "site", n_perm = 99, ambient = "flag")
  tf <- annotate_cmp_tf(fit, res, .toy_net(), n_perm = 100)
  tr <- cmp_triage(fit, res, context = "ctx", replicate_by = "site", tf = tf)
  expect_equal(nrow(tr), 1L)
  expect_equal(tr$n_genotype_markers, 0L)
  expect_true(tr$label %in% c("candidate", "candidate_associated"))
  expect_gte(tr$min_concordance, 0.99)
  expect_false(is.na(tr$context_p))
  expect_match(tr$tf_shared, "TFS\\+")
  nz <- c(genotype = "^G0[0-5][0-9]$", vdj = "^NONE$", sex = "^NONE$", ieg = "^NONE$")
  tr2 <- cmp_triage(fit, res, nuisance = nz)
  expect_equal(tr2$label, "genotype_candidate")
  expect_equal(tr2$n_genotype_markers, 12L)
  tr4 <- cmp_triage(fit, res, nuisance = c(platelet_rbc = "^G0[0-5][0-9]$"))
  expect_equal(tr4$label, "platelet_rbc")
  tr3 <- cmp_triage(fit, res, nuisance = nz, min_marker_genes = 13L)
  expect_match(tr3$flags, "genotype_marker_single")
})
