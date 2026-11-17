.toy <- function(n = 40, g = 60, seed = 3) {
  set.seed(seed)
  ids <- paste0("s", seq_len(n))
  batch <- rep(c("a", "b"), length.out = n)
  age <- rnorm(n)
  z <- rnorm(n)  # shared latent program, planted in both cell types
  mk <- function() {
    Y <- matrix(rnorm(n * g, sd = 0.5), n, g, dimnames = list(ids, paste0("g", seq_len(g))))
    Y[, 1:10] <- Y[, 1:10] + outer(z, rep(1, 10)) * 2
    Y + outer(ifelse(batch == "a", 2, -2), rnorm(g)) + outer(age, rnorm(g))
  }
  as_cell_program_data(list(A = mk(), B = mk()),
                       sample_metadata = data.frame(observation_id = ids, batch = batch, age = age))
}

# Repeated-visit toy: every cell type carries a patient effect shared across that
# patient's visits, but there is no cross-cell-type coupling (complete null).
.cmp_visit_toy <- function(n_pat = 40, visits = c("T1", "T2"), k = 3, seed = 1, rho = 0.8, cts = c("A", "B", "C")) {
  set.seed(seed)
  pat <- rep(sprintf("p%03d", seq_len(n_pat)), each = length(visits))
  vis <- rep(visits, n_pat)
  ids <- paste(pat, vis, sep = "_")
  n <- length(ids)
  scores <- loadings <- list()
  for (ct in cts) {
    U <- matrix(rnorm(n_pat * k), n_pat, k)
    S <- sqrt(rho) * U[match(pat, unique(pat)), , drop = FALSE] + sqrt(1 - rho) * matrix(rnorm(n * k), n, k)
    dimnames(S) <- list(ids, paste0(ct, "_", seq_len(k)))
    scores[[ct]] <- S
    loadings[[ct]] <- matrix(rnorm(50 * k), 50, k, dimnames = list(sprintf("G%02d", 1:50), colnames(S)))
  }
  structure(list(scores = scores, loadings = loadings, cell_types = cts, observation_ids = ids,
                 sample_metadata = data.frame(observation_id = ids, patient = pat, visit = vis,
                                              site = ifelse(match(pat, unique(pat)) %% 2 == 1, "s1", "s2"),
                                              visit_site = ifelse(seq_len(n) %% 4 < 2, "s1", "s2"),
                                              stringsAsFactors = FALSE)),
            class = c("cell_program_fit", "cell_program_data", "list"))
}

.pa_golden_fit <- function() {
  sim <- simulate_pseudobulk(n_samples = 40, n_genes = 200, program_size = 30, seed = 1)
  md <- data.frame(observation_id = paste0("obs", 1:40), site = rep(c("s1", "s2"), 20))
  x <- as_cell_program_data(sim$pseudobulk, sample_metadata = md)
  canonicalize_programs(fit_celltype_programs(x, max_factors = 3, seed = 1))
}
