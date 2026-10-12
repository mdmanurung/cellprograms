## registry.R -- named candidates (cmp-bench-v1). Each entry: function(fit, strata, seed) -> clusters df / list(clusters)
source("/mnt/shared-workspace/cmp/bench/candidates.R")
source("/mnt/shared-workspace/cmp/bench/cand_spec.R")
source("/mnt/shared-workspace/cmp/bench/cand_edge.R")
source("/mnt/shared-workspace/cmp/bench/cand_nodeagg.R")
source("/mnt/shared-workspace/cmp/bench/cand_ext.R")
source("/mnt/shared-workspace/cmp/bench/cand_flash.R")
source("/mnt/shared-workspace/cmp/bench/annot_lr.R")
source("/mnt/shared-workspace/cmp/bench/ambient.R")
NP <- 199L
REG <- list(
  ## F1 repo baselines
  b_pa_cc        = function(fit, strata, seed) run_method("pa_cc", fit, strata, NP, seed),
  b_node_cor     = function(fit, strata, seed) run_method("node_cor", fit, strata, NP, seed),
  b_subspace_sum = function(fit, strata, seed) run_method("subspace_sum", fit, strata, NP, seed),
  b_gene_match   = function(fit, strata, seed) run_method("gene_match", fit, strata, NP, seed),
  ## sanity
  s_shuffled     = function(fit, strata, seed) cand_edge(fit, strata, seed, stat = "spearman", shuffle = TRUE, n_perm = NP),
  ## F2/F3/F4 edge family
  e_pear_maxT_gr = function(fit, strata, seed) cand_edge(fit, strata, seed, stat = "pearson", n_perm = NP),
  e_rg_maxT_gr   = function(fit, strata, seed) cand_edge(fit, strata, seed, stat = "rankgauss", n_perm = NP),
  e_win_maxT_gr  = function(fit, strata, seed) cand_edge(fit, strata, seed, stat = "winsor", n_perm = NP),
  e_rg_bh_gr     = function(fit, strata, seed) cand_edge(fit, strata, seed, stat = "rankgauss", control = "bh", q = 0.05, n_perm = NP),
  e_rg_maxT_ld   = function(fit, strata, seed) cand_edge(fit, strata, seed, stat = "rankgauss", cluster = "leiden", n_perm = NP),
  e_rg_maxT_gr_nosite = function(fit, strata, seed) cand_edge(fit, strata, seed, stat = "rankgauss", resid_site = FALSE, n_perm = NP),
  ## F5 joint latent on node scores
  j_loo_mcca     = function(fit, strata, seed) cand_loo_mcca(fit, strata, seed, n_perm = NP),
  j_node_spec    = function(fit, strata, seed) cand_node_spec(fit, strata, seed, n_perm = NP),
  j_loo_mcca_rg  = function(fit, strata, seed) cand_loo_mcca(fit, strata, seed, n_perm = NP, rg = TRUE),
  j_node_spec_rg = function(fit, strata, seed) cand_node_spec(fit, strata, seed, n_perm = NP, rg = TRUE),
  j_flash        = function(fit, strata, seed) cand_flash_nodes(fit, strata, seed),
  e_rg_maxT_gr_amb3 = function(fit, strata, seed) cand_edge(fit, strata, seed, stat = "rankgauss", n_perm = NP, amb_thr = 0.3),
  e_rf_maxT_gr_amb3 = function(fit, strata, seed) cand_edge(fit, strata, seed, stat = "rankgauss", n_perm = NP, amb_thr = 0.3, rank_first = TRUE),
  e_rf_sd_gr_amb3 = function(fit, strata, seed) cand_edge(fit, strata, seed, stat = "rankgauss", control = "sdmaxT", n_perm = NP, amb_thr = 0.3, rank_first = TRUE),
  e_rf_ext = function(fit, strata, seed) cand_extend(fit, strata, seed, n_perm = NP),
  n_agg_ext = function(fit, strata, seed) cand_extend(fit, strata, seed, n_perm = NP, stage1 = cand_node_agg(fit, strata, seed, tau = 2, link_thr = 3, n_perm = NP)),
  n_agg_ms15 = function(fit, strata, seed) cand_node_agg(fit, strata, seed, tau = 2, link_thr = 3, n_perm = NP, min_shared = 15L),
  n_agg_ms30 = function(fit, strata, seed) cand_node_agg(fit, strata, seed, tau = 2, link_thr = 3, n_perm = NP, min_shared = 30L),
  n_agg_lt25 = function(fit, strata, seed) cand_node_agg(fit, strata, seed, tau = 2, link_thr = 2.5, n_perm = NP),
  n_agg_lt35 = function(fit, strata, seed) cand_node_agg(fit, strata, seed, tau = 2, link_thr = 3.5, n_perm = NP),
  n_agg_amb2 = function(fit, strata, seed) cand_node_agg(fit, strata, seed, tau = 2, link_thr = 3, n_perm = NP, amb_thr = 0.2),
  n_agg_amb5 = function(fit, strata, seed) cand_node_agg(fit, strata, seed, tau = 2, link_thr = 3, n_perm = NP, amb_thr = 0.5),
  n_agg_noanc = function(fit, strata, seed) cand_node_agg(fit, strata, seed, tau = 2, link_thr = 3, n_perm = NP, anchors = FALSE),
  n_agg_rg = function(fit, strata, seed) cand_node_agg(fit, strata, seed, tau = 2, link_thr = 3, n_perm = NP, rank_first = FALSE),
  n_agg_noamb = function(fit, strata, seed) cand_node_agg(fit, strata, seed, tau = 2, link_thr = 3, n_perm = NP, amb_thr = Inf),
  n_agg_minp = function(fit, strata, seed) cand_node_agg(fit, strata, seed, tau = 2, link_thr = 3, n_perm = NP, minp = TRUE),
  n_agg_minp_lt35 = function(fit, strata, seed) cand_node_agg(fit, strata, seed, tau = 2, link_thr = 3.5, n_perm = NP, minp = TRUE),
  n_agg_t2 = function(fit, strata, seed) cand_node_agg(fit, strata, seed, tau = 2, link_thr = 3, n_perm = NP),
  n_agg_t25 = function(fit, strata, seed) cand_node_agg(fit, strata, seed, tau = 2.5, link_thr = 3, n_perm = NP),
  e_rf_maxT_gr = function(fit, strata, seed) cand_edge(fit, strata, seed, stat = "rankgauss", n_perm = NP, rank_first = TRUE),
  e_rg_maxT_ld_amb3 = function(fit, strata, seed) cand_edge(fit, strata, seed, stat = "rankgauss", cluster = "leiden", n_perm = NP, amb_thr = 0.3),
  j_flash_perm   = function(fit, strata, seed) cand_flash_nodes(fit, strata, seed, perm = 39L)
)
