## Pipeline: stage-1 recovery, borrowing (L1/L2), integration, stability

test_that("stage-1 recovers Scenario 4 program structure", {
  r <- get_sim_fit(4L)
  expect_equal(vapply(r$fit$loadings, ncol, integer(1L)),
               c(B = 2L, CD4 = 1L, CD8 = 2L, NK = 2L, Mono = 2L))
  tz <- r$sim$truth$activities
  expect_gt(best_cor(r$fit, "B", tz$B_priv$z), 0.95)
  expect_gt(best_cor(r$fit, "Mono", tz$Mono_priv$z), 0.95)
  expect_gt(best_cor(r$fit, "CD8", tz$CD8NK_shared$z), 0.95)
  expect_gt(best_cor(r$fit, "NK", tz$CD8NK_shared$z), 0.95)
  expect_gt(best_cor(r$fit, "CD4", tz$IFN_global$z), 0.95)
})

test_that("scores-space sharing null is calibrated (Scenarios 1 and 7)", {
  s1 <- suppressWarnings(sharing_spectrum(get_sim_fit(1L)$fit, n_perm = 50L, seed = 1L))
  expect_true(all(s1$summary$n_shared_05 == 0L))
  s7 <- suppressWarnings(sharing_spectrum(get_sim_fit(7L)$fit, n_perm = 50L, seed = 1L))
  expect_true(all(s7$summary$n_shared_05 == 0L))
})

test_that("scores-space sharing spectrum recovers Scenario 4 structure", {
  ss <- suppressWarnings(sharing_spectrum(get_sim_fit(4L)$fit, n_perm = 50L, seed = 1L))
  summ <- ss$summary
  ## every pair shares IFN_global; CD8-NK additionally shares CD8NK_shared
  expect_true(all(summ$n_shared_05 >= 1L))
  expect_equal(summ$n_shared_05[summ$pair == "CD8 vs NK"], 2L)
})

test_that("L1 global prior pooling converges and does not degrade recovery", {
  r <- get_sim_fit(4L)
  f1 <- fit_celltype_programs(r$sim$x, max_factors = 10L, share_prior = "global",
                              seed = 42L)
  l1 <- f1$borrowing$level1
  expect_true(l1$n_iter <= 5L)
  expect_equal(vapply(f1$loadings, ncol, integer(1L)),
               c(B = 2L, CD4 = 1L, CD8 = 2L, NK = 2L, Mono = 2L))
  expect_gt(best_cor(f1, "CD4", r$sim$truth$activities$IFN_global$z), 0.95)
})

test_that("L2 mismatched-gene donor initialization shrinks to zero gracefully", {
  r <- get_sim_fit(4L)
  f2 <- fit_celltype_programs(r$sim$x, max_factors = 10L, seed = 42L,
                              init_from = list(NK = list(fit = r$fit, cell_type = "CD8")))
  rep <- f2$borrowing$level2$init_map$NK
  expect_equal(rep$n_borrowed, 2L)
  ## recovery unchanged: borrowed factors (wrong genes) must not distort
  expect_gt(best_cor(f2, "NK", r$sim$truth$activities$IFN_global$z), 0.95)
  expect_gt(best_cor(f2, "NK", r$sim$truth$activities$CD8NK_shared$z), 0.95)
})

test_that("PVA finds exactly the IFN-global and CD8-NK axes on Scenario 4", {
  r <- get_sim_fit(4L)
  f <- integrate_programs(r$fit, method = "pva", n_perm = 100L, seed = 1L)
  ax <- f$integration$axes
  expect_length(ax, 2L)
  S <- f$integration$states
  cors <- sapply(r$sim$truth$activities, function(a)
    apply(S, 2, function(s) abs(stats::cor(s, a$z, use = "pairwise.complete.obs"))))
  expect_gt(max(cors[, "IFN_global"]), 0.95)
  expect_gt(max(cors[, "CD8NK_shared"]), 0.8)
  ## traceability returns genes for the axis' cell types
  tr <- trace_state(f, state = "PV01", top_n = 5L)
  expect_true(length(tr) >= 2L)
})

test_that("block joint (scores level) classifies modules per ground truth", {
  r <- get_sim_fit(4L)
  blk <- fit_joint_block(r$fit, level = "scores", order = "shared_first", seed = 1L)
  cls <- blk$summary$class
  expect_true("global" %in% cls)
  expect_true("partial:CD8+NK" %in% cls)
  expect_equal(nrow(blk$summary), 2L)
})

test_that("EV-BIDIFAC (scores level) recovers global and partial modules", {
  r <- get_sim_fit(4L)
  evb <- fit_joint_evbidifac(r$fit, level = "scores")
  cls <- evb$summary$class
  expect_true("global" %in% cls)
  expect_true(any(grepl("CD8", cls) & grepl("NK", cls) & grepl("partial", cls)))
  S <- evb$states
  cors <- sapply(r$sim$truth$activities, function(a)
    apply(S, 2, function(s) abs(stats::cor(s, a$z, use = "pairwise.complete.obs"))))
  expect_gt(max(cors[, "IFN_global"], na.rm = TRUE), 0.95)
  expect_gt(max(cors[, "CD8NK_shared"], na.rm = TRUE), 0.95)
})

test_that("stability frequencies are in [0, 1] and programs are stable", {
  r <- get_sim_fit(4L)
  st <- assess_program_stability(r$fit, n_boot = 5L, seed = 7L, progress = FALSE)
  pp <- st$stability$per_program
  expect_true(all(pp$recovery_freq >= 0 & pp$recovery_freq <= 1))
  expect_true(all(pp$recovery_freq >= 0.8))
})

test_that("match_programs recovers identity on a same-data refit", {
  r <- get_sim_fit(4L)
  f_b <- canonicalize_programs(fit_celltype_programs(r$sim$x, max_factors = 10L,
                                                     seed = 99L))
  mm <- match_programs(r$fit, f_b, cell_type = "NK", method = "subspace",
                       sim_threshold = 0.5)
  expect_equal(nrow(mm$matches), 2L)
  expect_true(all(mm$matches$program_a == mm$matches$program_b))
  expect_true(all(abs(mm$matches$cosine) > 0.9))
})
