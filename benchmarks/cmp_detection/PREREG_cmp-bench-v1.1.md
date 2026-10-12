# PREREG cmp-bench-v1.1 -- held-out test of coordinated multicellular program (CMP) detection

Frozen: 2026-10-10, before any test-set candidate evaluation (test fits exist; no method has been run on them).

## Finalist (frozen)
`n_agg_noanc` = cand_node_agg(tau = 2, link_thr = 3, n_perm = 199, alpha = 0.05, min_shared = 20,
amb_thr = 0.3, rank_first = TRUE, anchors = FALSE, resid_site = TRUE): rank-gauss scores -> residualize on
strata -> drop ambient-like programs -> node statistic pooled over partner cell types -> max-node permutation
FWER (within-strata donor permutation) -> link significant nodes with |T| >= 3 -> greedy, one program per cell type.

## Comparators (frozen, reported, not selected among)
b_node_cor, b_pa_cc (repo), e_rf_maxT_gr_amb3 (edge max-T), n_agg_t2 (node+edge anchors), n_agg_minp_lt35
(fallback), b_gene_match (negative control), s_shuffled (sanity).

## Test set
Seeds 1001-1050 per base (param, semi=Stephenson) x 8 scenarios = 800 datasets (reduced from planned 100/base
for compute; decided before evaluation). Fits: main arm (point-Laplace, var_type 1, 2000 variable genes,
residualize on institute, K <= 10). Evaluator eval_cmp.R v1.1 (within mode permutes within site x context).

## Primary endpoint
membership_F1_pooled over positive scenarios (cmp_distinct, cmp_mixed, context, ambient, lr, weak), marginal mode.

## Success criteria for the finalist (all must hold)
1. F1 exceeds both repo baselines (b_node_cor, b_pa_cc): paired bootstrap 95% CI of the per-dataset difference > 0.
2. Null FWER (any call on a null dataset) <= 0.10 in each base (point estimate; Clopper-Pearson upper reported).
3. Ambient decoy rate <= 0.10 in each base.
4. Within-context false-call rate <= 0.10 in each base.
5. genes_only false-call rate <= 0.10 in each base.
6. s_shuffled F1 < 0.02 (harness sanity).
Failure of any criterion is reported as such; no re-tuning on the test set.

## Frozen file hashes (sha256)
996e8b6afd3cb401515c89f7da1af6a552f74df9e02803438a16f26c6102048c  common.R
16a5c8f73928e2c8f039edfff8963fdab580630a9b0b51096fcfd0fa46a72253  sim_cmp.R
ab940d3cc7c1b42c91a9993b0a913094cca8ef7627eaab2e7298c86fcfd4f9ff  config.R
49921921e5c860a8eaff705ad9abd418594c4c3fae03dd518d2a27b4126da36b  eval_cmp.R
cba72ad30db995a71c0363ca705e91ad3710df677de4df364df52be9f748daab  candidates.R
71e653d8fff9bde22239dfcb2de517972860307b697c15752fef51876d4c136a  cand_spec.R
e5aecb3abcb529b0cf22151a14b65c43e1aec0f5932733e27cbe0e64dcbc923c  cand_edge.R
5af2f7405d4d5e882ed4da890ffcd30f0aa1f0ae202ea146b74c09e01cd96eed  cand_nodeagg.R
7ace21a84237a70cce0b169d8c234d28c2be02a2e9a477b0d0cb94693080fa27  cand_ext.R
a84bad81405990603ee915258dc3770ba3714b0bb1b73904dd076fc66cdd85cd  cand_flash.R
ae6ff2ae27a532d8d965abd17f3daacac75a187461a14c4e437d543ad4231558  ambient.R
51dfd54f1169a734ac9ff65bff66cb8af638da27e2545da0eb2fdf098b6536f1  registry.R
8f65b8dfe09cbe9c03f63d7f9fb1efc6d017623d8cf5681766c8e18e5d5c90a7  run_cand.R
debf8293841a8d8bc0f71e6ba9ae4ef6eb726c4bb1a3e19b4c4506160493a2d5  fit_cache.R
