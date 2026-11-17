## probe_ebmf2_stage2.R -- does a second EBMF on the stacked stage-1 program scores find coordinated programs?
## Open question from the 2026-10-11 HANDOFF: on a toy, stage-2 EBMF found a 3-cell-type program in 10 of 12
## seeds but not a 2-cell-type program when 22% of donors were missing. This script re-asks it on the
## prereg-coord simulator and writes a TABLE of continuous recovery metrics. It makes no found / not-found
## call: that needs a threshold, which is the user's choice (impactful), and no figure has been approved.
## Result status: unsupported: no inspected figure. table fallback: the figure plan is not approved yet.
##
## stage 1: package defaults (fit_celltype_programs, canonicalize_programs).
## stage 2: flashier::flash on donors x nodes (all stage-1 programs, NA where a donor lacks the cell type),
##          package default prior (point_laplace) and var_type = 1, Kmax = STAGE2_KMAX.
## metric : for each planted program P and each stage-2 factor f,
##            share(P, f)   = sum |loading of P's best-matching node in each member ct| / sum |all loadings of f|
##            n_ct(P, f)    = number of P's member cell types whose best-matching node has |loading| >= 0.1 * max|loading of f|
##          per P the row with the largest share is kept. "Best-matching node" = node with the highest |cor|
##          to the planted activity. No threshold decides recovery.
## run: Rscript probe_ebmf2_stage2.R <out_csv> <seeds, comma-separated> [n_donors] [G] [amp]
a <- commandArgs(TRUE)
here <- local({ f <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]); dirname(normalizePath(f)) })
root <- normalizePath(file.path(here, "..", ".."))
out_csv <- a[1]; seeds <- as.integer(strsplit(a[2], ",")[[1]])
n_donors <- if (length(a) >= 3) as.integer(a[3]) else 100L
G <- if (length(a) >= 4) as.integer(a[4]) else 2000L
amp <- if (length(a) >= 5) as.numeric(a[5]) else stop("amp has no default; pass it as the fifth argument")
STAGE2_KMAX <- 10L
library(stats); library(utils)
source(file.path(here, "load_package.R")); load_package(root)
source(file.path(here, "sim_coord.R"))

one_seed <- function(seed) {
  d <- gen_coord("act_only", amp = amp, seed = seed, G = G, n_donors = n_donors, breadth = "all")
  fit <- canonicalize_programs(fit_celltype_programs(d$x, max_factors = 10, seed = seed))
  ids <- fit$observation_ids
  nodes <- do.call(cbind, lapply(fit$cell_types, function(ct) {
    S <- fit$scores[[ct]]
    if (!ncol(S)) return(NULL)
    M <- matrix(NA_real_, length(ids), ncol(S), dimnames = list(ids, colnames(S)))
    M[rownames(S), ] <- S
    M
  }))
  node_ct <- sub("_.*$", "", colnames(nodes))
  s2 <- flashier::flash(nodes, greedy_Kmax = STAGE2_KMAX, ebnm_fn = ebnm::ebnm_point_laplace,
                        var_type = 1, backfit = TRUE, verbose = 0)
  Lf <- s2$F_pm                                   # nodes x stage-2 factors
  rows <- list()
  for (P in d$truth$programs) {
    best <- vapply(P$members, function(ct) {       # best-matching node of this member cell type
      cols <- which(node_ct == ct); if (!length(cols)) return(NA_integer_)
      r <- abs(suppressWarnings(cor(nodes[names(P$z[[ct]]), cols, drop = FALSE], P$z[[ct]], use = "pairwise.complete.obs")))
      cols[which.max(r)]
    }, 1L)
    best <- best[!is.na(best)]
    if (!ncol(Lf) || !length(best)) { rows[[length(rows) + 1L]] <- data.frame(seed = seed, program = P$id, n_members = length(P$members),
      n_factors = ncol(Lf), share = NA_real_, n_ct = NA_integer_, stringsAsFactors = FALSE); next }
    sh <- apply(abs(Lf), 2, function(l) sum(l[best]) / max(sum(l), 1e-12))
    f <- which.max(sh)
    n_ct <- sum(abs(Lf[best, f]) >= 0.1 * max(abs(Lf[, f])))
    rows[[length(rows) + 1L]] <- data.frame(seed = seed, program = P$id, n_members = length(P$members),
      n_factors = ncol(Lf), share = unname(sh[f]), n_ct = n_ct, stringsAsFactors = FALSE)
  }
  do.call(rbind, rows)
}
res <- do.call(rbind, lapply(seeds, one_seed))
write.csv(res, out_csv, row.names = FALSE)
print(res, row.names = FALSE)
