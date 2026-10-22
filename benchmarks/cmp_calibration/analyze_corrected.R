# Bounded historical node-count summary. Does not certify calibration.
# Requested grid must enumerate every scenario/seed, including missing jobs.
summarize_node_counts <- function(results, requested) {
  need <- c("scenario", "seed")
  if (!all(need %in% names(requested)) || !nrow(requested) ||
      anyNA(requested[need]) || anyDuplicated(requested[need])) {
    stop("Requested grid needs unique nonmissing scenario/seed rows.")
  }
  cols <- c(need, "status", "null_n", "null_n_sig", "plant_n", "plant_n_sig")
  if (!all(cols %in% names(results)) || anyDuplicated(results[need]) ||
      anyNA(results[need])) stop("Results lack unique required node-count rows.")
  keys <- function(d) paste(d$scenario, d$seed, sep = "\r")
  if (length(setdiff(keys(results), keys(requested)))) stop("Unrequested result row.")
  r <- results[match(keys(requested), keys(results)), , drop = FALSE]
  r$scenario <- requested$scenario; r$seed <- requested$seed
  valid <- function(n, k) is.finite(n) & is.finite(k) & n >= 0 & k >= 0 &
    n == floor(n) & k == floor(k) & k <= n
  bad <- !is.na(r$status) & r$status == "ok" &
    (!valid(r$null_n, r$null_n_sig) | !valid(r$plant_n, r$plant_n_sig))
  if (any(bad, na.rm = TRUE)) stop("Invalid count in successful row.")
  upper <- function(k, n) if (!n) NA_real_ else if (k == n) 1 else
    stats::qbeta(.95, k + 1, n - k)
  out <- lapply(unique(requested$scenario), function(sc) {
    d <- r[r$scenario == sc, , drop = FALSE]
    # A result with no known null nodes cannot certify a null-family event.
    ok <- !is.na(d$status) & d$status == "ok" &
      valid(d$null_n, d$null_n_sig) & d$null_n > 0
    n <- nrow(d); nc <- sum(ok); nf <- sum(d$null_n_sig[ok] > 0)
    nn <- sum(d$null_n[ok]); fn <- sum(d$null_n_sig[ok])
    pn <- sum(d$plant_n[ok]); tp <- sum(d$plant_n_sig[ok])
    data.frame(scenario = sc, n_requested = n, n_callable = nc,
      n_failed_missing_or_uncallable = n - nc, n_any_false = nf,
      M06_conditional_node_family = if (nc) nf / nc else NA_real_,
      M06_conditional_upper95 = upper(nf, nc),
      M06_failure_inclusive_upper_rate = (nf + n - nc) / n,
      M06_failure_inclusive_upper95 = upper(nf + n - nc, n),
      M07_pooled_null_node_FPR = if (nn) fn / nn else NA_real_,
      M07_pooled_plant_power = if (pn) tp / pn else NA_real_,
      hypothesis_scope = "historical_known_score_node_counts",
      edge_or_set_FWER = NA_real_, calibration_certified = FALSE)
  })
  do.call(rbind, out)
}

if (sys.nframe() == 0L) {
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) != 3L) stop("Usage: Rscript analyze_corrected.R RESULTS.csv REQUESTED_GRID.csv NEW_OUTPUT.csv")
  if (file.exists(args[3L])) stop("Refuse to overwrite an existing/historical output.")
  out <- summarize_node_counts(utils::read.csv(args[1L]),
                                utils::read.csv(args[2L]))
  utils::write.csv(out, args[3L], row.names = FALSE)
}
