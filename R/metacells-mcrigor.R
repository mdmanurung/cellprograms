# mcRigor QC integration (native R; mcRigor is MIT, direct dependency).
#
# Liu & Li (2025) Nat Commun, doi:10.1038/s41467-025-63626-5.
# - mcDiv: divergence of the within-metacell feature-correlation matrix from
#   identity; null from a double permutation (within-cell preserving library
#   sizes, then within-feature); size-specific thresholds.
# - DubRate: fraction of dubious metacells; ZeroRate: fraction of zeros in the
#   metacell matrix; Score = 1 - w*DubRate - (1-w)*ZeroRate (w = 0.5 default).

.mcrigor_detect_group <- function(counts_g, membership, mode, mcrigor_args,
                                  group_id) {
    if (!requireNamespace("mcRigor", quietly = TRUE)) {
        cli::cli_abort(
            "{.pkg mcRigor} is required for metacell QC. Install it with {.code remotes::install_github('JSB-UCLA/mcRigor')}."
        )
    }
    # mcRigor's buildmc only aligns per-metacell metadata for NUMERIC
    # membership labels (its numeric branch strips Seurat's AggregateExpression
    # "g" prefix; string labels crash with "No cell overlap between new meta
    # data and Seurat object"). Map labels to integer codes here and remap the
    # normalized table back to the original metacell names below.
    lab <- as.character(membership)
    codes <- as.integer(factor(lab, levels = sort(unique(lab))))
    names(codes) <- names(membership)
    # minimal Seurat object holding the group's single-cell counts
    obj_g <- Seurat::CreateSeuratObject(
        counts = counts_g, assay = "RNA", min.features = 0, min.cells = 0
    )
    membership_df <- data.frame(
        matrix(codes, ncol = 1, dimnames = list(names(codes),
                                                     as.character(1))),
        stringsAsFactors = FALSE, check.names = FALSE
    )
    args <- list(
        Nrep = 50, gene_filter = 0.1, feature_use = 2000,
        test_cutoff = 0.01, cor_method = "pearson",
        aggregate_method = if (identical(mode, "average")) "mean" else "sum",
        draw = FALSE, step_save = FALSE, output_file = tempfile(fileext = ".rds"),
        # skip mcRigor's per-metacell covariate transfer: it crashes with
        # string metacell labels (Seurat name mangling misaligns value names
        # vs the metacell object cells) and is not needed for the QC
        # statistics (mcDiv / DubRate / ZeroRate / Score)
        fields = character(0)
    )
    args <- utils::modifyList(args, mcrigor_args)
    base_args <- list(
        obj_singlecell = obj_g,
        cell_membership = membership_df,
        tgamma = as.character(1)
    )
    res <- do.call(mcRigor::mcRigor_DETECT, utils::modifyList(base_args, args))
    out <- .mcrigor_normalize(res, group_id = group_id)
    # remap integer codes back to the original metacell labels
    code_to_lab <- sort(unique(lab))
    code_idx <- suppressWarnings(as.integer(out$metacell))
    if (!any(is.na(code_idx)) && max(code_idx) <= length(code_to_lab)) {
        out$metacell <- code_to_lab[code_idx]
    }
    out
}

.mcrigor_normalize <- function(res, group_id) {
    # Normalization of mcRigor_DETECT output into a per-metacell data.frame
    # with columns: metacell, mcDiv, mcDiv_threshold, dubious.
    # Pinned against mcRigor (Liu & Li 2025) return structure:
    #   res$TabMC: gamma, ZeroRate, size, T_org (observed divergence),
    #     T_colperm, T_rowperm*, T_bothperm*, TT_div, depctr_TT_div,
    #     mcRigor (per-metacell verdict string, e.g. "trustworthy")
    #   res$thre: data.frame(size, thre) with metacell ids as rownames
    #     (size-specific divergence threshold)
    #   res$mc_res: named chr verdict per metacell
    tab <- NULL
    if (is.data.frame(res)) {
        tab <- res
    } else if (is.list(res)) {
        for (nm in c("TabMC", "tab", "result", "results", "mcDiv")) {
            if (!is.null(res[[nm]]) && is.data.frame(res[[nm]])) {
                tab <- res[[nm]]
                break
            }
        }
        if (is.null(tab)) {
            dfs <- Filter(is.data.frame, res)
            if (length(dfs) > 0) tab <- dfs[[1]]
        }
    }
    if (is.null(tab)) {
        stop("unrecognized mcRigor_DETECT return structure")
    }
    tab <- as.data.frame(tab)
    if ("metacell" %in% colnames(tab)) {
        mc_col <- tab$metacell
    } else {
        mc_col <- rownames(tab)
    }
    pick <- function(cands) {
        hit <- intersect(cands, colnames(tab))
        if (length(hit) > 0) tab[[hit[1]]] else rep(NA_real_, nrow(tab))
    }
    out <- data.frame(
        metacell = as.character(mc_col),
        Group = group_id,
        mcDiv = as.numeric(pick(c("mcDiv", "T_org", "divergence", "score"))),
        mcDiv_threshold = as.numeric(pick(c("threshold", "cutoff", "mcDiv_threshold"))),
        dubious = as.logical(pick(c("dubious", "is_dubious", "flag"))),
        stringsAsFactors = FALSE
    )
    # size-specific threshold from res$thre (metacell ids as rownames)
    if (all(is.na(out$mcDiv_threshold)) && !is.null(res$thre) &&
            is.data.frame(res$thre) && "thre" %in% colnames(res$thre)) {
        idx <- match(out$metacell, rownames(res$thre))
        out$mcDiv_threshold <- as.numeric(res$thre$thre[idx])
    }
    # per-metacell verdict string ("trustworthy" / otherwise dubious)
    if (all(is.na(out$dubious))) {
        verdict <- pick(c("mcRigor", "verdict", "mc_res"))
        if (is.character(verdict) && !all(is.na(verdict))) {
            # NA verdicts stay NA: `!NA %in% ...` would coerce them to TRUE
            out$dubious <- ifelse(is.na(verdict), NA,
                                  !verdict %in% c("trustworthy", "TRUE", "good"))
        }
    }
    if (all(is.na(out$dubious)) && !all(is.na(out$mcDiv)) &&
            !all(is.na(out$mcDiv_threshold))) {
        out$dubious <- out$mcDiv > out$mcDiv_threshold
    }
    out
}

#' Optimize SEACells granularity with mcRigor
#'
#' Runs SEACells at several granularities (cells per metacell) per group and
#' calls `mcRigor::mcRigor_OPTIMIZE` to select the granularity maximizing the
#' mcRigor Score (1 - w*DubRate - (1-w)*ZeroRate). Expensive: one SEACells fit
#' per group per candidate granularity.
#'
#' @param seurat_obj A Seurat object.
#' @param group.by Metadata column(s) defining groups (as in
#'   [MetacellsByGroupsSEACells()]).
#' @param gammas Candidate granularities as cells-per-metacell values.
#'   Default `c(25, 50, 100)` (mcRigor semi-synthetic optimum was 50 for
#'   SEACells).
#' @param reduction, dims, assay, slot, layer, n_top_genes, n_pc,
#'   n_waypoint_eigs, waypt_proportion, n_neighbors, convergence_epsilon,
#'   max_iter, use_sparse, min_cells, seed, verbose
#'   Passed to [MetacellsByGroupsSEACells()] internals.
#' @param mcrigor_args Optional named list overriding mcRigor_OPTIMIZE
#'   arguments (`D_bw`, `optim_method`, `weight`, `Nrep`, ...).
#'
#' @return A list with `optimal_gamma` (named list, per group), `scores` (per
#'   group x gamma Score table), and `membership` (named list of per-group
#'   membership vectors for the chosen granularity, one entry per cell; NULL
#'   when no group produced a result). Gammas whose membership collapses to a
#'   single metacell cannot be scored by mcRigor and are dropped with a
#'   warning (`NA` in `scores`).
#'
#' @export
OptimizeGammaSEACells <- function(
    seurat_obj,
    group.by = "seurat_clusters",
    gammas = c(25, 50, 100),
    reduction = "pca",
    dims = NULL,
    assay = "RNA",
    slot = "counts",
    layer = "counts",
    n_top_genes = 1500,
    n_pc = 50,
    n_waypoint_eigs = 10,
    waypt_proportion = 0.9,
    n_neighbors = 15,
    convergence_epsilon = 1e-5,
    max_iter = 20,
    use_sparse = TRUE,
    min_cells = 100,
    mcrigor_args = list(),
    seed = 20261016,
    verbose = TRUE
) {
    if (!requireNamespace("mcRigor", quietly = TRUE)) {
        cli::cli_abort(
            "{.pkg mcRigor} is required for granularity optimization. Install it with {.code remotes::install_github('JSB-UCLA/mcRigor')}."
        )
    }
    if (!requireNamespace("Seurat", quietly = TRUE)) {
        cli::cli_abort(
            "{.pkg Seurat} is required for granularity optimization. Install it with {.code install.packages('Seurat')}."
        )
    }
    group_ids <- .mc_get_group_ids(seurat_obj, group.by)
    levels_g <- sort(unique(group_ids))

    membership_cols <- list()
    score_rows <- list()
    optimal <- list()

    opt_root <- tempfile(pattern = "cpmc_opt_")   # unique per run (no same-second collisions)
    for (gi in seq_along(levels_g)) {
        g <- levels_g[gi]
        cells_g <- names(group_ids)[group_ids == g]
        if (length(cells_g) < min_cells) {
            warning("group '", g, "' skipped: ", length(cells_g),
                    " cells < min_cells = ", min_cells, call. = FALSE)
            next
        }
        if (verbose) message("[OptimizeGamma] group '", g, "': ",
                             length(gammas), " granularities")
        counts_g <- .get_counts_matrix(seurat_obj, assay, slot, layer, cells_g)
        emb_g <- if (!is.null(reduction)) {
            .get_reduction_embeddings(seurat_obj, reduction, dims, cells_g)
        } else {
            NULL
        }
        mem_cols <- list()
        for (k in seq_along(gammas)) {
            gamma <- gammas[k]
            n_seacells <- .mc_n_seacells_for_group(length(cells_g), gamma, NULL)
            out_dir <- file.path(opt_root, paste0("g", gi, "_", k))
            paths <- .export_group_for_python(counts_g, emb_g, out_dir)
            params <- list(
                n_top_genes = n_top_genes, n_pc = n_pc,
                n_waypoint_eigs = n_waypoint_eigs,
                waypt_proportion = waypt_proportion, n_neighbors = n_neighbors,
                convergence_epsilon = convergence_epsilon, max_iter = max_iter,
                use_sparse = use_sparse, mode = "sum"
            )
            .run_seacells_in_basilisk(paths, out_dir, n_seacells, params, seed)
            mc_out <- .read_metacell_outputs(out_dir)
            lab <- .hard_labels_vector(mc_out, cells_g,
                                       paste0("mc_", seq_len(ncol(mc_out$mtx))))
            mem_cols[[k]] <- lab
        }
        membership_df <- as.data.frame(
            mem_cols, stringsAsFactors = FALSE, check.names = FALSE
        )
        colnames(membership_df) <- as.character(gammas)
        rownames(membership_df) <- cells_g

        # mcRigor's buildmc aligns per-metacell sizes through Seurat's
        # AggregateExpression colnames, which mangle membership labels
        # (digit-leading values gain a "g" prefix; "_" becomes "-" under
        # Seurat v5) and only undo the digit case (an is.numeric gate) —
        # the "mc_<k>" labels used here would never realign, failing every
        # group with "No cell overlap between new meta data and Seurat
        # object". Score on integer codes per gamma column (the same
        # contract as .mcrigor_detect_group) and keep the labels for the
        # return value. A gamma whose membership collapsed to a single
        # metacell cannot be scored (no between-metacell variance, and
        # AggregateExpression ignores 1-value groupings) — drop it up
        # front with a warning instead of failing the whole group.
        coded <- .mc_optimize_membership_codes(mem_cols, gammas, cells_g)
        if (length(coded$drop_gamma)) {
            warning("group '", g, "': gamma ", paste(coded$drop_gamma, collapse = ", "),
                    " yield(s) a single metacell; mcRigor cannot score it — dropped",
                    call. = FALSE)
        }
        if (!length(coded$keep_gamma)) {
            warning("group '", g, "': no candidate granularity produced >= 2 metacells",
                    call. = FALSE)
            next
        }

        obj_g <- Seurat::CreateSeuratObject(
            counts = counts_g, assay = "RNA", min.features = 0, min.cells = 0
        )
        args <- list(Nrep = 50, draw = FALSE, step_save = FALSE,
                     output_file = tempfile(fileext = ".rds"))
        args <- utils::modifyList(args, mcrigor_args)
        res <- tryCatch(
            do.call(
                mcRigor::mcRigor_OPTIMIZE,
                utils::modifyList(
                    list(obj_singlecell = obj_g, cell_membership = coded$code_df),
                    args
                )
            ),
            error = function(e) {
                warning("mcRigor_OPTIMIZE failed for group '", g, "': ",
                        conditionMessage(e))
                NULL
            }
        )
        if (is.null(res)) next
        # normalize: expect per-gamma scores; pin exact structure at smoke test
        best <- .mcrigor_optimize_pick(res, coded$keep_gamma)
        optimal[[g]] <- best$gamma
        # keep the full gamma grid in the scores table; dropped gammas score NA
        sc <- rep(NA_real_, length(gammas))
        sc[match(as.character(coded$keep_gamma), as.character(gammas))] <- best$scores
        score_rows[[g]] <- data.frame(group = g, gamma = gammas, score = sc)
        membership_cols[[g]] <- membership_df[[as.character(best$gamma)]]
    }

    list(
        optimal_gamma = optimal,
        scores = do.call(rbind, score_rows),
        membership = if (length(membership_cols) > 0) membership_cols else NULL
    )
}

.mc_optimize_membership_codes <- function(mem_cols, gammas, cells_g) {
    # Prepare OptimizeGammaSEACells' membership table for mcRigor_OPTIMIZE:
    # drop gammas whose membership collapsed to a single metacell (mcRigor
    # cannot score them) and integer-code the remaining columns so buildmc's
    # size alignment survives Seurat's AggregateExpression label mangling
    # ("g" prefix on digit-leading values, "_" -> "-" under Seurat v5).
    n_unique <- vapply(mem_cols, function(l) length(unique(l)), integer(1))
    keep_idx <- which(n_unique >= 2L)
    keep_gamma <- gammas[keep_idx]
    drop_gamma <- gammas[n_unique < 2L]
    code_df <- if (length(keep_idx)) {
        df <- as.data.frame(
            lapply(mem_cols[keep_idx],
                   function(l) as.integer(factor(l, levels = unique(l)))),
            stringsAsFactors = FALSE, check.names = FALSE
        )
        colnames(df) <- as.character(keep_gamma)
        rownames(df) <- cells_g
        df
    } else {
        NULL
    }
    list(keep_gamma = keep_gamma, drop_gamma = drop_gamma, code_df = code_df)
}

.mcrigor_optimize_pick <- function(res, gammas) {
    # Defensive extraction of per-gamma scores from mcRigor_OPTIMIZE output.
    tab <- if (is.data.frame(res)) res else {
        dfs <- Filter(is.data.frame, res)
        if (length(dfs) > 0) dfs[[1]] else NULL
    }
    if (is.null(tab)) {
        stop("unrecognized mcRigor_OPTIMIZE return structure")
    }
    score_col <- intersect(c("Score", "score", "opt_score"), colnames(tab))[1]
    gamma_col <- intersect(c("gamma", "Gamma", "granularity"), colnames(tab))[1]
    if (is.na(score_col)) {
        # fall back: assume rows follow gammas order
        scores <- as.numeric(tab[[1]])
    } else {
        scores <- as.numeric(tab[[score_col]])
    }
    if (!is.na(gamma_col)) {
        gvals <- tab[[gamma_col]]
        aligned <- scores[match(as.character(gammas), as.character(gvals))]
        # Keep the alignment whenever at least one gamma matched; the previous
        # any(NA) fallback discarded correct matches and restored raw order,
        # misaligning every score after a partial match.
        if (!all(is.na(aligned))) scores <- aligned
    }
    if (length(scores) != length(gammas)) {
        scores <- rep(NA_real_, length(gammas))
    }
    best_idx <- if (all(is.na(scores))) NA_integer_ else which.max(scores)
    list(gamma = gammas[best_idx], scores = scores)
}

#' Access metacell results stored by MetacellsByGroupsSEACells
#'
#' @param seurat_obj A Seurat object processed by
#'   [MetacellsByGroupsSEACells()].
#' @param wgcna_name Name under which results were stored. Default
#'   `"SEACells"`.
#'
#' @return The stored list: `metacell_object`, `soft_assignments`, `qc`,
#'   `fit_stats`, `zerorate_overall`, `dubrate_overall`, `score_overall`,
#'   `params`, `skipped_groups`, `call`.
#'
#' @export
GetMetacellSEACells <- function(seurat_obj, wgcna_name = "SEACells") {
    out <- seurat_obj@misc[[wgcna_name]]
    if (is.null(out)) {
        stop("no metacell results under misc[['", wgcna_name,
             "']]; run MetacellsByGroupsSEACells first")
    }
    out
}
