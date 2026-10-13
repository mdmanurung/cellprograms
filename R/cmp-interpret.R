# Interpretation helpers for coordinated multicellular programs: TF footprints
# of member programs (annotate_cmp_tf) and a one-row-per-CMP triage table
# (cmp_triage). Neither affects whether a CMP is called.

#' Format an OmniPath CollecTRI table as a signed TF-target network
#'
#' Follows the decoupler convention: one row per TF-target pair, mode of
#' regulation (`mor`) -1 when the consensus is inhibition only and +1
#' otherwise (stimulation, both, or unknown); NF-kB and AP-1 complex members
#' are merged into `NFKB` and `AP1`.
#'
#' @param raw Data frame from the OmniPath `interactions` endpoint with
#'   `datasets=collectri&genesymbols=yes` (columns `source`,
#'   `source_genesymbol`, `target_genesymbol`, `consensus_stimulation`,
#'   `consensus_inhibition`).
#' @return Data frame with `source`, `target`, `mor`.
#' @examples
#' # `raw` mimics the OmniPath interactions endpoint (datasets=collectri)
#' raw <- data.frame(
#'   source = c("P1", "COMPLEX:P2_P3", "P4"),
#'   source_genesymbol = c("STAT1", "JUN_FOS", "MYC"),
#'   target_genesymbol = c("IRF1", "MMP1", "CDKN1A"),
#'   consensus_stimulation = c(TRUE, TRUE, FALSE),
#'   consensus_inhibition = c(FALSE, FALSE, TRUE)
#' )
#' prepare_collectri(raw)
#' @export
prepare_collectri <- function(raw) {
  need <- c("source", "source_genesymbol", "target_genesymbol", "consensus_stimulation", "consensus_inhibition")
  miss <- setdiff(need, names(raw))
  if (length(miss)) .stopf("`raw` lacks column(s): %s.", paste(miss, collapse = ", "))
  tf <- as.character(raw$source_genesymbol)
  cx <- grepl("^COMPLEX:", raw$source)
  parts <- strsplit(tf, "_")
  nfkb <- c("NFKB1", "NFKB2", "REL", "RELA", "RELB")
  ap1 <- c("JUN", "JUNB", "JUND", "FOS", "FOSB", "FOSL1", "FOSL2")
  tf[cx & vapply(parts, function(p) all(p %in% nfkb), TRUE)] <- "NFKB"
  tf[cx & vapply(parts, function(p) all(p %in% ap1), TRUE)] <- "AP1"
  .lgl <- function(x) as.logical(toupper(as.character(x)))
  inh <- .lgl(raw$consensus_inhibition) %in% TRUE & !(.lgl(raw$consensus_stimulation) %in% TRUE)
  out <- data.frame(source = tf, target = as.character(raw$target_genesymbol), mor = ifelse(inh, -1, 1),
                    stringsAsFactors = FALSE)
  out <- out[!duplicated(out[, c("source", "target")]), , drop = FALSE]
  rownames(out) <- NULL
  out
}

# signed, CMP-oriented loading of one member program (optionally shrunk by 1 - lfsr)
.cmp_oriented_loading <- function(fit, ct, prog, sg = 1, use_lfsr = TRUE) {
  W <- fit$loadings[[ct]]
  if (is.null(W) || !prog %in% colnames(W)) return(NULL)
  w <- W[, prog]
  lf <- fit$fits[[ct]]$flash$F_lfsr
  if (use_lfsr && !is.null(lf) && nrow(lf) == length(w)) w <- w * (1 - lf[, match(prog, colnames(W))])
  if (is.finite(sg)) w <- sg * w
  w
}

.cmp_std_cols <- function(X) {
  X <- sweep(X, 2, colMeans(X))
  s <- sqrt(colSums(X^2))
  sweep(X, 2, ifelse(s > 0, s, 1), "/")
}

#' Transcription-factor footprints of coordinated programs
#'
#' For each CMP member program, a univariate linear model (decoupler-style
#' ULM) of its signed loadings on each TF's signed regulon (`mor` for targets,
#' 0 otherwise) within the program's gene universe. Loadings are multiplied by
#' the member's CMP sign, so a positive `t` means the TF's targets go up in
#' donors with a high CMP score; with `use_lfsr`, loadings are shrunk by
#' (1 - lfsr). `p` is the ULM t-test p-value (what a gene-label permutation
#' converges to, without its 1/(n_perm + 1) floor, which would otherwise stop
#' BH from reaching significance over hundreds of TFs), `padj` is BH within
#' the member, and `p_fwer` uses the gene-permutation maximum of |t| over all
#' TFs tested for that member. This competitive test asks whether targets are
#' more (or less) loaded than random genes. Because loadings have a natural
#' zero, a regulon of unloaded genes in a program skewed to one sign would
#' look "depleted"; a TF is therefore called (`sig`) only if it also passes a
#' self-contained sign-flip test (`p_self < q_self`: the mor-signed target
#' loadings are coherently non-zero) in the same direction as `t`.
#'
#' Members are then summarized per CMP and TF. A TF is `"shared"` if it is
#' called (`sig`) with the same sign in at least two members and at
#' least half of the members it was tested in; `"specific"` if significant in
#' exactly one member; and `"mixed"` otherwise (significant in several members
#' with discordant signs or in fewer than half). `z` is the Stouffer
#' combination of member t-statistics (converted to normal scores), used for
#' ranking only: members share genes and priors, so it is not a calibrated
#' test.
#'
#' The gene-permutation null treats genes as exchangeable. Co-regulated genes
#' violate that assumption, so `padj` is anticonservative when a large block of
#' correlated genes drives a program; large, promiscuous regulons (e.g. NFKB,
#' AP1, SP1, MYC) reach significance more readily. A footprint says which
#' regulons a program resembles, not which TF caused the coordination.
#'
#' @param fit A canonicalized `cell_program_fit`.
#' @param cmp A `cmp_result` or its `cmps` table.
#' @param net Data frame with `source` (TF), `target` and `mor` (+1/-1), e.g.
#'   from `prepare_collectri()`.
#' @param min_size Minimum regulon targets within the program's genes (default 5).
#' @param n_perm Permutations for `p_fwer` and `p_self` (default 1000; 0 for normal approximations).
#' @param q Threshold on the competitive `padj` (default 0.05).
#' @param q_self Threshold on the self-contained `p_self` (default 0.05).
#' @param use_lfsr Shrink loadings by (1 - lfsr) when available (default TRUE).
#' @param seed Random seed.
#' @return List of class `cmp_tf`: `members` (one row per member and tested
#'   TF: `cmp`, `cell_type`, `program_id`, `tf`, `n_targets`, `t`, `p`,
#'   `padj`, `p_fwer`, `p_self`, `mean_target_loading` (mean mor-signed target
#'   loading, scaled by the member's maximum |loading|), `sig`, `tf_loading` =
#'   the TF gene's own oriented loading scaled
#'   by the member's maximum |loading|, NA if absent) and `summary` (one row per
#'   CMP and TF significant in at least one member: `cmp`, `tf`, `n_tested`,
#'   `n_up`, `n_down`, `z`, `class`, `members_sig`).
#' @examples
#' # Toy canonicalized fit: program 1 of cell types A and B shares donor
#' # activity z (a planted coordinated program) and top genes G1-G15.
#' set.seed(1)
#' ids <- paste0("d", seq_len(40))
#' genes <- paste0("G", seq_len(100))
#' z <- rnorm(40)
#' mk <- function(ct) {
#'   S <- matrix(rnorm(40 * 2), 40, 2, dimnames = list(ids, paste0(ct, "_", 1:2)))
#'   S[, 1] <- z + rnorm(40, sd = 0.3)
#'   W <- matrix(rnorm(length(genes) * 2, sd = 0.1), length(genes), 2,
#'               dimnames = list(genes, colnames(S)))
#'   W[1:15, 1] <- W[1:15, 1] + 1
#'   list(S = S, W = W)
#' }
#' ta <- mk("A"); tb <- mk("B")
#' fit <- structure(list(
#'   scores = list(A = ta$S, B = tb$S), loadings = list(A = ta$W, B = tb$W),
#'   cell_types = c("A", "B"), observation_ids = ids,
#'   sample_metadata = data.frame(observation_id = ids,
#'                                group = rep(c("ctrl", "stim"), each = 20))),
#'   class = c("cell_program_fit", "cell_program_data", "list"))
#' res <- find_coordinated_programs(fit, n_perm = 99, ambient = "none")
#' net <- rbind(data.frame(source = "TF1", target = paste0("G", 1:15), mor = 1),
#'              data.frame(source = "TF2", target = paste0("G", 20:34), mor = 1))
#' tf <- annotate_cmp_tf(fit, res, net, n_perm = 99)
#' tf$summary
#' @export
annotate_cmp_tf <- function(fit, cmp, net, min_size = 5L, n_perm = 1000L, q = 0.05, q_self = 0.05,
                            use_lfsr = TRUE, seed = 1L) {
  tab <- .cmp_table(cmp)
  if (!all(c("source", "target", "mor") %in% names(net))) .stopf("`net` needs columns source, target, mor.")
  if (!nrow(tab)) return(structure(list(members = NULL, summary = NULL), class = "cmp_tf"))
  set.seed(seed)
  tfs_all <- unique(net$source)
  mem <- lapply(seq_len(nrow(tab)), function(i) {
    ct <- tab$cell_type[i]; pr <- tab$program_id[i]
    sg <- if (!is.null(tab$sign) && is.finite(tab$sign[i])) tab$sign[i] else 1
    w <- .cmp_oriented_loading(fit, ct, pr, sg, use_lfsr)
    if (is.null(w)) return(NULL)
    g <- names(w)
    nn <- net[net$target %in% g, , drop = FALSE]
    sz <- table(nn$source)
    keep <- names(sz)[sz >= min_size]
    if (!length(keep)) return(NULL)
    nn <- nn[nn$source %in% keep, , drop = FALSE]
    M <- matrix(0, length(g), length(keep), dimnames = list(g, keep))
    M[cbind(match(nn$target, g), match(nn$source, keep))] <- nn$mor
    M <- M[, apply(M, 2, function(x) length(unique(x)) > 1L), drop = FALSE]
    if (!ncol(M)) return(NULL)
    n <- length(w); df <- n - 2
    Ms <- .cmp_std_cols(M)
    tstat <- function(R) R * sqrt(df / pmax(1 - R^2, 1e-12))
    t <- tstat(drop(crossprod(.cmp_std_cols(matrix(w)), Ms)))
    s_self <- drop(crossprod(w, M))                               # sum over targets of mor * loading
    p <- 2 * stats::pt(-abs(t), df)   # ULM p; a gene-permutation p converges to this but is floored at 1/(n_perm+1)
    if (n_perm > 0) {
      P <- vapply(seq_len(n_perm), function(b) sample(w), numeric(n))
      Tp <- abs(tstat(crossprod(.cmp_std_cols(P), Ms)))          # n_perm x TF
      mx <- apply(Tp, 1, max)
      p_fwer <- vapply(abs(t), function(a) (1 + sum(mx >= a)) / (n_perm + 1), 1)
      Fl <- matrix(sample(c(-1, 1), n * n_perm, replace = TRUE), n, n_perm) * w
      Sp <- abs(crossprod(Fl, M))                                  # sign-flip null of s_self
      p_self <- (1 + colSums(sweep(Sp, 2, abs(s_self), ">="))) / (n_perm + 1)
    } else {
      p_fwer <- pmin(1, 2 * stats::pt(-abs(t), df) * length(t))
      p_self <- 2 * stats::pnorm(-abs(s_self) / sqrt(pmax(drop(crossprod(w^2, M^2)), 1e-24)))
    }
    tfl <- w[match(colnames(M), g)] / max(abs(w), 1e-12)
    data.frame(cmp = tab$cmp[i], cell_type = ct, program_id = pr, tf = colnames(M),
               n_targets = colSums(M != 0), t = unname(t), p = unname(p),
               padj = stats::p.adjust(p, "BH"), p_fwer = unname(p_fwer), p_self = unname(p_self),
               mean_target_loading = unname(s_self) / colSums(M != 0) / max(abs(w), 1e-12),
               tf_loading = unname(tfl),
               z = -sign(unname(t)) * stats::qnorm(stats::pt(-abs(unname(t)), df, log.p = TRUE), log.p = TRUE),
               stringsAsFactors = FALSE)
  })
  mem <- do.call(rbind, mem)
  if (is.null(mem)) return(structure(list(members = NULL, summary = NULL), class = "cmp_tf"))
  rownames(mem) <- NULL
  mem$sig <- mem$padj < q & mem$p_self < q_self & sign(mem$t) == sign(mem$mean_target_loading)
  sig <- mem[mem$sig, , drop = FALSE]
  keys <- unique(sig[, c("cmp", "tf")])
  summ <- lapply(seq_len(nrow(keys)), function(j) {
    a <- mem[mem$cmp == keys$cmp[j] & mem$tf == keys$tf[j], , drop = FALSE]
    s <- a[a$sig, , drop = FALSE]
    up <- sum(s$t > 0); dn <- sum(s$t < 0); nt <- nrow(a)
    cls <- if (max(up, dn) >= 2L && min(up, dn) == 0L && max(up, dn) >= nt / 2) "shared" else
      if (nrow(s) == 1L) "specific" else "mixed"
    data.frame(cmp = keys$cmp[j], tf = keys$tf[j], n_tested = nt, n_up = up, n_down = dn,
               z = sum(a$z) / sqrt(nt), class = cls,
               members_sig = paste0(s$cell_type, ifelse(s$t > 0, "+", "-"), collapse = ","),
               stringsAsFactors = FALSE)
  })
  summ <- do.call(rbind, summ)
  if (!is.null(summ)) {
    summ <- summ[order(summ$cmp, match(summ$class, c("shared", "mixed", "specific")), -abs(summ$z)), , drop = FALSE]
    rownames(summ) <- NULL
  }
  mem$z <- NULL
  structure(list(members = mem, summary = summ), class = "cmp_tf")
}

#' Print a CMP transcription-factor triage result
#'
#' @param x A `cmp_tf` from [annotate_cmp_tf()].
#' @param ... Ignored.
#' @return `x`, invisibly.
#' @export
print.cmp_tf <- function(x, ...) {
  s <- x$summary
  cat(sprintf("<cmp_tf> %d member-TF tests; %d CMP-TF footprints (%d shared, %d specific, %d mixed)\n",
              if (is.null(x$members)) 0L else nrow(x$members), if (is.null(s)) 0L else nrow(s),
              sum(s$class == "shared"), sum(s$class == "specific"), sum(s$class == "mixed")))
  invisible(x)
}

#' Default nuisance-gene patterns for CMP triage
#'
#' Regular expressions for gene classes that produce cross-cell-type donor
#' co-variation without coordinated biology: immunoglobulin/TCR V(D)J and
#' constant genes (`vdj`), genes with common germline deletion, copy-number or
#' HLA-haplotype expression variation (`genotype`), sex-chromosome genes
#' (`sex`), immediate-early/dissociation stress genes (`ieg`) and platelet or
#' erythrocyte contamination markers (`platelet_rbc`; these cells are usually
#' not fitted, so the ambient score cannot see them). The
#' `genotype` list is a screening heuristic, not validated against genotype
#' data; extend or replace it for your cohort.
#'
#' @return Named character vector of regular expressions.
#' @examples
#' cmp_nuisance_genes()
#' grep(cmp_nuisance_genes()["vdj"], c("IGHV1-2", "IGKC", "GAPDH"), value = TRUE)
#' @export
cmp_nuisance_genes <- function() {
  c(vdj = "^(IG[HKL]V|IG[KL]C|IGH[GAMDE]|JCHAIN|TR[ABGD]V)",
    genotype = "^(GSTM1|GSTT1|GSTT2B?|UTS2|ERAP2|CCL3L1|CCL3L3|CCL4L2|HLA-DRB[345]|HLA-DQ[AB][12])$",
    sex = "^(XIST|TSIX|RPS4Y1|DDX3Y|UTY|KDM5D|EIF1AY|USP9Y|ZFY|TXLNGY|NLGN4Y)$",
    ieg = "^(FOS|FOSB|JUN|JUNB|EGR1|EGR2|EGR3|ATF3|IER2|IER3|NR4A1|DUSP1|PPP1R15A|KLF6)$",
    platelet_rbc = "^(PF4|PPBP|GP9|GP1BB|TUBB1|ITGA2B|GNG11|CAVIN2|TREML1|PTCRA|HBB|HBA1|HBA2|HBD|ALAS2)$")
}

#' One-row-per-CMP triage table
#'
#' Collects the checks needed before interpreting a CMP: the recurrent top
#' genes of its members and how many of them fall in nuisance classes
#' (`cmp_nuisance_genes()`), the fraction of ambient-flagged members, an
#' optional association with a donor variable, optional sign replication of
#' member-pair co-variation within levels of a grouping variable (e.g. site),
#' and optional TF footprints from `annotate_cmp_tf()`.
#'
#' Recurrent top genes: the `top_n` genes with the largest CMP-oriented
#' loadings in each member, ranked by how many members share them; the first
#' `n_recurrent` are kept. Flags (all thresholds are arguments):
#' `ambient` (ambient fraction >= `max_ambient_frac`), `platelet_rbc`
#' (>= `min_marker_genes` contamination markers), `vdj` (V(D)J fraction of
#' recurrent genes >= `max_vdj_frac`), `sex_chromosome` and
#' `genotype_candidate` (>= `min_marker_genes` marker genes), `stress_ieg`
#' (>= 3 IEGs), plus soft notes `genotype_marker_single` (one marker gene) and
#' `low_replication` (minimum within-level concordance < `min_concordance`).
#' `label` is the first hard flag in that order, else `"candidate"`, or
#' `"candidate_associated"` when the context association has p < 0.05. Labels
#' are heuristics for prioritization, not classifications.
#'
#' @param fit A canonicalized `cell_program_fit`.
#' @param cmp A `cmp_result` or its `cmps` table (needs `ambient_flag` for the
#'   ambient check; use `ambient = "flag"` in `find_coordinated_programs()`).
#' @param context Optional donor variable (column name or named vector) for
#'   `cmp_context_association()`.
#' @param covariates Optional `sample_metadata` columns to adjust for.
#' @param replicate_by Optional grouping variable (column name or named vector),
#'   e.g. site. Within each level with >= `min_level_donors` donors, member-pair
#'   Pearson correlations of CMP-oriented scores are computed; concordance is
#'   the fraction that are positive.
#' @param tf Optional `cmp_tf` result.
#' @param nuisance Named regex vector (default `cmp_nuisance_genes()`).
#' @param top_n,n_recurrent Top genes per member (15) and recurrent genes kept (12).
#' @param max_ambient_frac,max_vdj_frac,min_marker_genes,min_concordance Flag thresholds.
#' @param min_level_donors Minimum donors per `replicate_by` level (default 15).
#' @return Data frame, one row per CMP.
#' @examples
#' # Toy canonicalized fit: program 1 of cell types A and B shares donor
#' # activity z (a planted coordinated program) and top genes G1-G15.
#' set.seed(1)
#' ids <- paste0("d", seq_len(40))
#' genes <- paste0("G", seq_len(100))
#' z <- rnorm(40)
#' mk <- function(ct) {
#'   S <- matrix(rnorm(40 * 2), 40, 2, dimnames = list(ids, paste0(ct, "_", 1:2)))
#'   S[, 1] <- z + rnorm(40, sd = 0.3)
#'   W <- matrix(rnorm(length(genes) * 2, sd = 0.1), length(genes), 2,
#'               dimnames = list(genes, colnames(S)))
#'   W[1:15, 1] <- W[1:15, 1] + 1
#'   list(S = S, W = W)
#' }
#' ta <- mk("A"); tb <- mk("B")
#' fit <- structure(list(
#'   scores = list(A = ta$S, B = tb$S), loadings = list(A = ta$W, B = tb$W),
#'   cell_types = c("A", "B"), observation_ids = ids,
#'   sample_metadata = data.frame(observation_id = ids,
#'                                group = rep(c("ctrl", "stim"), each = 20))),
#'   class = c("cell_program_fit", "cell_program_data", "list"))
#' res <- find_coordinated_programs(fit, n_perm = 99, ambient = "none")
#' cmp_triage(fit, res, context = "group", replicate_by = "group")
#' @export
cmp_triage <- function(fit, cmp, context = NULL, covariates = NULL, replicate_by = NULL, tf = NULL,
                       nuisance = cmp_nuisance_genes(), top_n = 15L, n_recurrent = 12L,
                       max_ambient_frac = 0.5, max_vdj_frac = 0.5, min_marker_genes = 2L,
                       min_concordance = 0.75, min_level_donors = 15L) {
  tab <- .cmp_table(cmp)
  if (!nrow(tab)) return(data.frame())
  cnt <- function(g, cl) if (cl %in% names(nuisance)) g[grepl(nuisance[[cl]], g)] else character()
  rows <- lapply(unique(tab$cmp), function(k) {
    m <- tab[tab$cmp == k, , drop = FALSE]
    tops <- lapply(seq_len(nrow(m)), function(i) {
      sg <- if (!is.null(m$sign) && is.finite(m$sign[i])) m$sign[i] else 1
      w <- fit$loadings[[m$cell_type[i]]][, m$program_id[i]] * sg
      names(sort(-w))[seq_len(min(top_n, length(w)))]
    })
    tb <- sort(table(unlist(tops)), decreasing = TRUE)
    rec <- utils::head(names(tb), n_recurrent)
    gt <- cnt(rec, "genotype"); sx <- cnt(rec, "sex")
    data.frame(cmp = k, n_members = nrow(m), cell_types = paste(m$cell_type, collapse = ","),
               members = paste0(m$program_id, ifelse(!is.null(m$sign) & m$sign < 0, "(-)", ""), collapse = ","),
               recurrent_top_genes = paste(rec, collapse = ","),
               ambient_frac = if (is.null(m$ambient_flag)) NA_real_ else mean(m$ambient_flag %in% TRUE),
               vdj_frac = length(cnt(rec, "vdj")) / max(length(rec), 1L),
               n_genotype_markers = length(gt), genotype_markers = paste(gt, collapse = ","),
               n_sex_markers = length(sx), sex_markers = paste(sx, collapse = ","),
               n_ieg = length(cnt(rec, "ieg")), n_platelet_rbc = length(cnt(rec, "platelet_rbc")),
               stringsAsFactors = FALSE)
  })
  out <- do.call(rbind, rows)

  out$context_estimate <- out$context_p <- NA_real_
  if (!is.null(context)) {
    ca <- cmp_context_association(fit, tab, context, covariates)
    best <- ca[order(ca$p), , drop = FALSE]
    best <- best[!duplicated(best$cmp), , drop = FALSE]
    i <- match(as.character(out$cmp), as.character(best$cmp))
    out$context_p <- best$p[i]; out$context_estimate <- best$estimate[i]
  }

  out$min_concordance <- NA_real_; out$concordance_by_level <- NA_character_
  if (!is.null(replicate_by)) {
    rb <- .cmp_resolve_var(fit, replicate_by, "replicate_by", allow_na = TRUE)
    levs <- sort(unique(stats::na.omit(as.character(rb))))
    for (r in seq_len(nrow(out))) {
      m <- tab[tab$cmp == out$cmp[r], , drop = FALSE]
      if (nrow(m) < 2L) next
      V <- sapply(seq_len(nrow(m)), function(i) {
        v <- fit$scores[[m$cell_type[i]]][, m$program_id[i]]
        sg <- if (!is.null(m$sign) && is.finite(m$sign[i])) m$sign[i] else 1
        sg * v[match(fit$observation_ids, names(v))]
      })
      rownames(V) <- fit$observation_ids
      cc <- vapply(levs, function(l) {
        ids <- names(rb)[rb %in% l]
        X <- V[ids, , drop = FALSE]
        if (sum(stats::complete.cases(X)) < min_level_donors) return(NA_real_)
        C <- suppressWarnings(stats::cor(X, use = "pairwise.complete.obs"))
        mean(C[upper.tri(C)] > 0, na.rm = TRUE)
      }, 1)
      out$concordance_by_level[r] <- paste(sprintf("%s:%.2f", levs, cc), collapse = ";")
      out$min_concordance[r] <- suppressWarnings(min(cc, na.rm = TRUE))
    }
    out$min_concordance[!is.finite(out$min_concordance)] <- NA
  }

  out$tf_shared <- out$tf_specific <- NA_character_
  if (!is.null(tf) && !is.null(tf$summary)) {
    for (r in seq_len(nrow(out))) {
      s <- tf$summary[tf$summary$cmp == out$cmp[r], , drop = FALSE]
      sh <- s[s$class == "shared", , drop = FALSE]
      sh <- utils::head(sh[order(-abs(sh$z)), , drop = FALSE], 5)
      out$tf_shared[r] <- paste(sprintf("%s%s(%d/%d)", sh$tf, ifelse(sh$z > 0, "+", "-"),
                                        pmax(sh$n_up, sh$n_down), sh$n_tested), collapse = ",")
      sp <- s[s$class == "specific", , drop = FALSE]
      sp <- utils::head(sp[order(-abs(sp$z)), , drop = FALSE], 5)
      out$tf_specific[r] <- paste(sprintf("%s:%s", sp$tf, sp$members_sig), collapse = ",")
    }
  }

  hard <- cbind(ambient = !is.na(out$ambient_frac) & out$ambient_frac >= max_ambient_frac,
                platelet_rbc = out$n_platelet_rbc >= min_marker_genes,
                vdj = out$vdj_frac >= max_vdj_frac,
                sex_chromosome = out$n_sex_markers >= min_marker_genes,
                genotype_candidate = out$n_genotype_markers >= min_marker_genes,
                stress_ieg = out$n_ieg >= 3L)
  soft <- cbind(genotype_marker_single = out$n_genotype_markers > 0 & out$n_genotype_markers < min_marker_genes,
                low_replication = !is.na(out$min_concordance) & out$min_concordance < min_concordance)
  fl <- cbind(hard, soft)
  out$flags <- apply(fl, 1, function(x) paste(colnames(fl)[x], collapse = ","))
  out$label <- apply(hard, 1, function(x) if (any(x)) colnames(hard)[which(x)[1]] else NA_character_)
  cand <- is.na(out$label)
  out$label[cand] <- ifelse(!is.na(out$context_p[cand]) & out$context_p[cand] < 0.05,
                            "candidate_associated", "candidate")
  rownames(out) <- NULL
  out
}
