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
