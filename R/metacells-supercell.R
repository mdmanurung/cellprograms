# SuperCell 2.0 metacell backend (pure R; the benchmark-default path).
#
# Benchmark (benchmarks/metacell_comparison/, ledger 2026-10-22): SuperCell 2.0
# + mcRigor matches SEACells on dubious-detection ARI (gap 0.019 <= 0.05) and
# program-recovery fidelity (gap 5.4% <= 10%) at ~200x lower membership
# construction cost, so it is the default R-native metacell path. SEACells
# remains available unchanged (R/metacells-seacells.R).
#
# SuperCell (Bilous et al. 2022, doi:10.1093/bioinformatics/btac104; v2.0,
# GfellerLab/SuperCell@89b34c0, supercell-2.0 branch) is used via
# SCimplify_for_Seurat_v5: gamma is natively "cells per metacell", matching
# the SEACells interface's cells_per_metacell parameter.

.sc_gamma_for_target <- function(obj_g, target, reduction, dims, assay, k.knn,
                                 tol = 0.10, max_iter = 14) {
    # Bisect log-gamma to land within +/-10% of the target metacell count.
    # SCimplify gamma ~= cells per metacell (N_mc = floor(n_cells / gamma));
    # higher gamma -> fewer metacells. Benchmark-proven pattern (adaptive
    # bounds: lo = max(2, n/(target*4)), hi = n/(target/4)).
    n_cells <- ncol(obj_g)
    lo <- max(2, n_cells / (target * 4))
    hi <- n_cells / max(1, target / 4)
    best <- NULL
    for (it in seq_len(max_iter)) {
        g <- sqrt(lo * hi)
        sc <- try(
            SuperCell::SCimplify_for_Seurat_v5(
                obj_g, gamma = g, k.knn = k.knn, assay = assay,
                reduction = list(reduction), dims = list(dims),
                return.seurat = FALSE, verbose = FALSE
            ),
            silent = TRUE
        )
        if (inherits(sc, "try-error") || is.null(sc$membership)) {
            hi <- g
            next
        }
        n_mc <- length(unique(sc$membership))
        best <- list(gamma = g, n = n_mc, mem = sc$membership, iter = it)
        if (abs(n_mc - target) <= tol * target) break
        if (n_mc > target) lo <- g else hi <- g
    }
    best
}

.sc_run_simplify <- function(obj_g, gamma, reduction, dims, assay, k.knn) {
    sc <- SuperCell::SCimplify_for_Seurat_v5(
        obj_g, gamma = gamma, k.knn = k.knn, assay = assay,
        reduction = list(reduction), dims = list(dims),
        return.seurat = FALSE, verbose = FALSE
    )
    if (is.null(sc$membership) || !length(sc$membership)) {
        stop("SCimplify_for_Seurat_v5 returned no membership")
    }
    sc
}

.sc_aggregate_counts <- function(counts_g, mem_codes, mc_names, mode) {
    # Aggregate single-cell counts into metacells by sparse indicator product
    # (benchmark-proven: counts (genes x cells) %*% t(ind) (cells x k) ->
    # genes x k). mode = "sum" (raw sums, SEACells default) or "average"
    # (divide by metacell size).
    uid <- sort(unique(mem_codes))
    ind <- Matrix::sparseMatrix(
        i = match(mem_codes, uid), j = seq_along(mem_codes),
        x = rep(1, length(mem_codes)),
        dims = c(length(uid), length(mem_codes))
    )
    mc_mtx <- counts_g %*% Matrix::t(ind)   # sparse genes x metacells
    if (identical(mode, "average")) {
        sizes <- as.numeric(table(mem_codes)[uid])
        mc_mtx <- sweep(mc_mtx, 2, sizes, "/")
    }
    dimnames(mc_mtx) <- list(rownames(counts_g), mc_names)
    mc_mtx
}

.sc_prepare_reduction <- function(obj_g, reduction, dims, n_top_genes, n_pc,
                                  seed) {
    # Return list(reduction_name, dims_vector) for the group's object. When
    # reduction = NULL, compute a standard Seurat PCA on the group (the
    # benchmark's preprocessing: normalize -> HVGs -> scale -> PCA).
    if (!is.null(reduction)) {
        if (!(reduction %in% names(obj_g@reductions))) {
            stop("reduction '", reduction, "' not found; available: ",
                 paste(names(obj_g@reductions), collapse = ", "))
        }
        emb <- obj_g@reductions[[reduction]]@cell.embeddings
        d <- if (!is.null(dims)) dims else seq_len(ncol(emb))
        return(list(reduction = reduction, dims = d))
    }
    obj_g <- Seurat::NormalizeData(obj_g, verbose = FALSE)
    obj_g <- Seurat::FindVariableFeatures(obj_g, nfeatures = n_top_genes,
                                          verbose = FALSE)
    obj_g <- Seurat::ScaleData(obj_g, verbose = FALSE)
    obj_g <- Seurat::RunPCA(obj_g, npcs = n_pc, seed = seed, verbose = FALSE)
    list(reduction = "PCA", dims = seq_len(min(n_pc, ncol(obj_g@reductions[["PCA"]]@cell.embeddings))))
}

#' Construct metacells by groups with SuperCell and mcRigor QC
#'
#' Constructs transcriptomically homogeneous metacells from a Seurat object
#' using SuperCell 2.0 (Bilous et al. 2022) as the aggregation backend,
#' computed separately within each level of the supplied groupings (typically
#' donor, cell type, or condition), following the hdWGCNA `MetacellsByGroups`
#' interface (Morabito et al. 2023). This is the package's **default**,
#' pure-R metacell path: a preregistered benchmark
#' (`benchmarks/metacell_comparison/`) found SuperCell + mcRigor comparable to
#' SEACells on dubious-detection ARI and program-recovery fidelity at ~200x
#' lower construction cost. Metacell quality is governed by mcRigor (Liu & Li
#' 2025): per-metacell mcDiv dubious detection with double-permutation nulls,
#' DubRate/ZeroRate scoring, and optional granularity optimization via
#' [OptimizeGammaSuperCell()]. The SEACells backend remains available in
#' [MetacellsByGroupsSEACells()].
#'
#' @param seurat_obj A Seurat object with raw counts and a dimensionality
#'   reduction (default `pca`). Set `reduction = NULL` to compute a standard
#'   Seurat PCA per group on the fly.
#' @param group.by Character vector of metadata column names defining the
#'   stratification (e.g. `c("donor", "celltype")`). Metacells are constructed
#'   separately within each combination of levels.
#' @param ident.group Optional metadata column name used as the metacell-name
#'   prefix; defaults to the interaction of `group.by`.
#' @param reduction Name of the dimensionality reduction SCimplify builds the
#'   cell-to-cell graph on. Default `"pca"`; `NULL` computes a PCA per group.
#' @param dims Dimensions of `reduction` to use (indices); default all
#'   computed dimensions of the reduction.
#' @param assay Assay to aggregate. Default `RNA`.
#' @param slot Slot (Seurat v4) to extract counts from. Default `counts`.
#' @param layer Layer (Seurat v5) to extract counts from. Default `counts`.
#' @param mode How to aggregate constituent cells into metacells: `"sum"`
#'   (raw-count sums, the SEACells default and the benchmark configuration) or
#'   `"average"` (divide by metacell size). Default `"sum"`.
#' @param cells_per_metacell Target number of single cells per metacell  - 
#'   SuperCell's native `gamma`. Default 75 (SEACells guidance; hdWGCNA k
#'   guidance: 20-75).
#' @param target_metacells Optional cap on the number of metacells per group.
#'   When set, gamma is bisected (log scale) to land within +/-10% of the
#'   target; otherwise `gamma = cells_per_metacell` is used directly.
#' @param min_cells Groups with fewer than this many cells are skipped with a
#'   warning. Default 100.
#' @param k.knn Number of neighbors for the cell-to-cell KNN graph. Default
#'   30 (SCimplify default).
#' @param n_top_genes Number of highly variable genes for the on-the-fly PCA
#'   when `reduction = NULL`. Default 2000.
#' @param n_pc Number of PCs for the on-the-fly PCA when `reduction = NULL`.
#'   Default 30.
#' @param run_mcrigor If TRUE (default), run mcRigor dubious-metacell
#'   detection per group and attach mcDiv/DubRate/ZeroRate/Score to the QC
#'   table.
#' @param mcrigor_args Optional named list overriding mcRigor_DETECT
#'   arguments: `Nrep` (50), `gene_filter` (0.1), `feature_use` (2000),
#'   `test_cutoff` (0.01), `cor_method` ("pearson"), `aggregate_method`
#'   (defaults to `"sum"`/`"mean"` following `mode`).
#' @param seed Random seed (set before each group's SCimplify call; the graph
#'   construction is deterministic but the parameter is kept for API parity
#'   and on-the-fly PCA reproducibility). Default 20261016.
#' @param name Reserved for interface parity with
#'   [MetacellsByGroupsSEACells()] (metacell names are `<group>_<i>`).
#'   Default `"SC"`.
#' @param wgcna_name Name under which results are stored in
#'   `seurat_obj@misc`. Default `"SuperCell"`.
#' @param verbose Print progress. Default TRUE.
#'
#' @return The Seurat object with results stored in `seurat_obj@misc[[wgcna_name]]`:
#'   `metacell_object` (Seurat object of metacells, counts = aggregated
#'   expression, metadata = Group, n_cells, the individual `group.by` columns,
#'   and mcRigor QC columns), `soft_assignments` (always NULL for this
#'   backend: SCimplify returns hard memberships only), `qc` (per-metacell QC
#'   data.frame), `fit_stats` (per-group gamma and metacell count),
#'   `params` (all call parameters), `skipped_groups` (data.frame), and
#'   `call`. Metacells are named `<group>_<i>` and column-bound across groups.
#'
#' @details Groups are processed independently; each group's cells are
#'   subsetted (with their reduction), SCimplify partitions the KNN graph,
#'   metacell counts are summarized from raw counts, and mcRigor QC is
#'   computed against the group's single-cell counts. The result is
#'   hdWGCNA-interoperable in the same way as the SEACells backend.
#'
#' @references
#' Bilous et al. (2022) Bioinformatics, doi:10.1093/bioinformatics/btac104;
#' Liu & Li (2025) Nat Commun, doi:10.1038/s41467-025-63626-5;
#' Morabito et al. (2023) Cell Rep Methods, doi:10.1016/j.crmeth.2023.100498.
#'
#' @export
MetacellsByGroupsSuperCell <- function(
    seurat_obj,
    group.by = "seurat_clusters",
    ident.group = NULL,
    reduction = "pca",
    dims = NULL,
    assay = "RNA",
    slot = "counts",
    layer = "counts",
    mode = c("sum", "average"),
    cells_per_metacell = 75,
    target_metacells = NULL,
    min_cells = 100,
    k.knn = 30,
    n_top_genes = 2000,
    n_pc = 30,
    run_mcrigor = TRUE,
    mcrigor_args = list(),
    seed = 20261016,
    name = "SC",
    wgcna_name = "SuperCell",
    verbose = TRUE
) {
    mode <- match.arg(mode)
    if (!requireNamespace("Seurat", quietly = TRUE) ||
        !requireNamespace("SeuratObject", quietly = TRUE)) {
        cli::cli_abort(
            "{.pkg Seurat} and {.pkg SeuratObject} are required for metacell construction. Install them with {.code install.packages('Seurat')}."
        )
    }
    if (!requireNamespace("SuperCell", quietly = TRUE)) {
        cli::cli_abort(
            "{.pkg SuperCell} is required for the SuperCell metacell backend. Install it with {.code remotes::install_github('GfellerLab/SuperCell@89b34c078dba0b91289a5a23cf75a6db3a46d232')}."
        )
    }
    if (!methods::is(seurat_obj, "Seurat")) {
        stop("seurat_obj must be a Seurat object")
    }
    missing_cols <- setdiff(group.by, colnames(seurat_obj@meta.data))
    if (length(missing_cols) > 0) {
        stop("group.by columns not found in meta.data: ",
             paste(missing_cols, collapse = ", "))
    }
    clashing <- intersect(group.by, c("Group", "n_cells"))
    if (length(clashing) > 0) {
        cli::cli_abort("{.arg group.by} column(s) {.val {clashing}} clash with reserved metacell metadata names; rename the metadata column first.")
    }
    if (!is.null(reduction) && !(reduction %in% names(seurat_obj@reductions))) {
        stop("reduction '", reduction, "' not found; available: ",
             paste(names(seurat_obj@reductions), collapse = ", "))
    }
    if (isTRUE(run_mcrigor) && !requireNamespace("mcRigor", quietly = TRUE)) {
        cli::cli_warn("{.pkg mcRigor} is not installed; skipping dubious-metacell QC (QC columns will be NA). Install it with {.code remotes::install_github('JSB-UCLA/mcRigor')}.")
        run_mcrigor <- FALSE
    }

    group_ids <- .mc_get_group_ids(seurat_obj, group.by)
    levels_g <- sort(unique(group_ids))
    prefix_base <- if (!is.null(ident.group)) ident.group else paste(group.by, collapse = "_")

    params <- list(
        group.by = group.by, ident.group = ident.group, reduction = reduction,
        dims = dims, assay = assay, slot = slot, layer = layer, mode = mode,
        cells_per_metacell = cells_per_metacell,
        target_metacells = target_metacells, min_cells = min_cells,
        k.knn = k.knn, n_top_genes = n_top_genes, n_pc = n_pc,
        run_mcrigor = run_mcrigor, seed = seed, name = name
    )

    mc_list <- list(); obs_list <- list(); qc_list <- list()
    skipped <- list(); fit_stats_list <- list()

    for (gi in seq_along(levels_g)) {
        g <- levels_g[gi]
        cells_g <- names(group_ids)[group_ids == g]
        g_prefix <- .mc_sanitize(paste0(prefix_base, "_", g))
        if (verbose) message("[", gi, "/", length(levels_g), "] group '", g,
                             "': ", length(cells_g), " cells")
        if (length(cells_g) < min_cells) {
            warning("group '", g, "' skipped: ", length(cells_g),
                    " cells < min_cells = ", min_cells, call. = FALSE)
            skipped[[g]] <- data.frame(group = g, n_cells = length(cells_g),
                                       reason = "below min_cells")
            next
        }

        counts_g <- .get_counts_matrix(seurat_obj, assay, slot, layer, cells_g)
        obj_g <- subset(seurat_obj, cells = cells_g)
        red <- .sc_prepare_reduction(obj_g, reduction, dims, n_top_genes, n_pc,
                                     seed)

        # gamma selection: direct or bisected to the target
        gamma_used <- cells_per_metacell
        bisect_iter <- NA_integer_
        if (!is.null(target_metacells)) {
            hit <- .sc_gamma_for_target(obj_g, target_metacells, red$reduction,
                                        red$dims, assay, k.knn)
            if (is.null(hit) || is.null(hit$mem)) {
                warning("group '", g, "': gamma bisection found no membership ",
                        "for target_metacells = ", target_metacells,
                        "; skipped", call. = FALSE)
                skipped[[g]] <- data.frame(group = g, n_cells = length(cells_g),
                                           reason = "gamma bisection failed")
                next
            }
            gamma_used <- hit$gamma
            bisect_iter <- hit$iter
            mem_codes <- hit$mem
        } else {
            sc <- .sc_run_simplify(obj_g, gamma_used, red$reduction, red$dims,
                                   assay, k.knn)
            mem_codes <- sc$membership
        }
        mem_codes <- mem_codes[cells_g]   # deterministic cell order
        uid <- sort(unique(mem_codes))
        m <- length(uid)
        mc_names <- paste0(g_prefix, "_", seq_len(m))

        mc_mtx <- .sc_aggregate_counts(counts_g, mem_codes, mc_names, mode)
        mc_obs <- data.frame(
            metacell = mc_names, n_cells = as.integer(table(mem_codes)[uid]),
            stringsAsFactors = FALSE, row.names = mc_names
        )
        mc_obs$Group <- g
        for (col in group.by) {
            mc_obs[[col]] <- unique(as.character(seurat_obj@meta.data[cells_g, col]))[1]
        }

        # ---- mcRigor QC (native R, shared wrapper) ----
        qc_g <- NULL
        if (run_mcrigor) {
            mem_lab <- mc_names[match(mem_codes, uid)]
            names(mem_lab) <- cells_g
            qc_g <- tryCatch(
                .mcrigor_detect_group(
                    counts_g = counts_g,
                    membership = mem_lab,
                    mode = mode, mcrigor_args = mcrigor_args, group_id = g
                ),
                error = function(e) {
                    warning("mcRigor QC failed for group '", g, "': ",
                            conditionMessage(e))
                    NULL
                }
            )
        }
        if (is.null(qc_g)) {
            qc_g <- data.frame(
                metacell = mc_names, Group = g, n_cells = mc_obs$n_cells,
                mcDiv = NA_real_, mcDiv_threshold = NA_real_,
                dubious = NA, stringsAsFactors = FALSE
            )
        }
        if (!"n_cells" %in% colnames(qc_g)) {
            qc_g$n_cells <- mc_obs$n_cells[match(qc_g$metacell, mc_names)]
        }
        for (cc in setdiff(c("metacell", "Group", "n_cells", "mcDiv",
                             "mcDiv_threshold", "dubious"), colnames(qc_g))) {
            qc_g[[cc]] <- NA
        }
        zr_mc <- 1 - Matrix::colSums(mc_mtx != 0) / nrow(mc_mtx)
        qc_g$zerorate_mc <- as.numeric(zr_mc)[match(qc_g$metacell, mc_names)]
        qc_list[[g]] <- qc_g

        mc_list[[g]] <- mc_mtx
        obs_list[[g]] <- mc_obs
        fit_stats_list[[g]] <- data.frame(
            group = g, gamma = gamma_used, n_metacells = m,
            bisection_iterations = bisect_iter, stringsAsFactors = FALSE
        )
    }

    if (length(mc_list) == 0) {
        stop("no groups passed min_cells = ", min_cells)
    }

    mc_all <- Reduce(Matrix::cbind2, mc_list)
    obs_all <- do.call(rbind, obs_list)
    qc_all <- do.call(rbind, qc_list)
    fit_stats_all <- do.call(rbind, fit_stats_list)

    zerorate_overall <- 1 - Matrix::nnzero(mc_all) / length(mc_all)
    dubrate_overall <- if (any(!is.na(qc_all$dubious))) {
        mean(qc_all$dubious[!is.na(qc_all$dubious)] == TRUE)
    } else {
        NA_real_
    }
    score_overall <- if (!is.na(dubrate_overall)) {
        1 - 0.5 * dubrate_overall - 0.5 * zerorate_overall
    } else {
        NA_real_
    }

    m_obj <- Seurat::CreateSeuratObject(
        counts = mc_all, assay = assay, min.features = 0, min.cells = 0
    )
    m_obj[["Group"]] <- factor(obs_all$Group)
    m_obj[["n_cells"]] <- obs_all$n_cells
    for (col in group.by) {
        m_obj[[col]] <- obs_all[[col]]
    }
    qc_idx <- match(colnames(m_obj), qc_all$metacell)
    for (col in setdiff(colnames(qc_all), c("metacell", "Group", "n_cells", group.by))) {
        m_obj[[col]] <- qc_all[[col]][qc_idx]
    }

    skipped_df <- if (length(skipped) > 0) do.call(rbind, skipped) else
        data.frame(group = character(0), n_cells = integer(0), reason = character(0))

    out <- list(
        metacell_object = m_obj,
        soft_assignments = NULL,
        qc = qc_all,
        fit_stats = fit_stats_all,
        zerorate_overall = zerorate_overall,
        dubrate_overall = dubrate_overall,
        score_overall = score_overall,
        params = params,
        skipped_groups = skipped_df,
        call = match.call()
    )
    seurat_obj@misc[[wgcna_name]] <- out
    invisible(seurat_obj)
}

#' Optimize SuperCell granularity with mcRigor
#'
#' Runs SuperCell's SCimplify at several granularities (cells per metacell  - 
#' SuperCell's native `gamma`) per group and calls `mcRigor::mcRigor_OPTIMIZE`
#' to select the granularity maximizing the mcRigor Score
#' (1 - w*DubRate - (1-w)*ZeroRate). Cheap relative to the SEACells analogue:
#' one graph partition per group per candidate granularity (seconds).
#'
#' @param seurat_obj A Seurat object.
#' @param group.by Metadata column(s) defining groups (as in
#'   [MetacellsByGroupsSuperCell()]).
#' @param gammas Candidate granularities as cells-per-metacell values.
#'   Default `c(25, 50, 100)`.
#' @param reduction, dims, assay, slot, layer, n_top_genes, n_pc, k.knn,
#'   min_cells, seed, verbose
#'   Passed to [MetacellsByGroupsSuperCell()] internals.
#' @param mcrigor_args Optional named list overriding mcRigor_OPTIMIZE
#'   arguments (`Nrep`, `weight`, ...).
#'
#' @return A list with `optimal_gamma` (named list, per group), `scores` (per
#'   group x gamma Score table), and `membership` (named list of per-group
#'   membership vectors for the chosen granularity, one entry per cell; NULL
#'   when no group produced a result). Gammas whose membership collapses to a
#'   single metacell cannot be scored by mcRigor and are dropped with a
#'   warning (`NA` in `scores`).
#'
#' @export
OptimizeGammaSuperCell <- function(
    seurat_obj,
    group.by = "seurat_clusters",
    gammas = c(25, 50, 100),
    reduction = "pca",
    dims = NULL,
    assay = "RNA",
    slot = "counts",
    layer = "counts",
    min_cells = 100,
    k.knn = 30,
    n_top_genes = 2000,
    n_pc = 30,
    mcrigor_args = list(),
    seed = 20261016,
    verbose = TRUE
) {
    if (!requireNamespace("mcRigor", quietly = TRUE)) {
        cli::cli_abort(
            "{.pkg mcRigor} is required for granularity optimization. Install it with {.code remotes::install_github('JSB-UCLA/mcRigor')}."
        )
    }
    if (!requireNamespace("SuperCell", quietly = TRUE)) {
        cli::cli_abort(
            "{.pkg SuperCell} is required for the SuperCell metacell backend. Install it with {.code remotes::install_github('GfellerLab/SuperCell@89b34c078dba0b91289a5a23cf75a6db3a46d232')}."
        )
    }
    group_ids <- .mc_get_group_ids(seurat_obj, group.by)
    levels_g <- sort(unique(group_ids))
    gammas <- .mc_validate_gammas(gammas)

    membership_cols <- list()
    score_rows <- list()
    optimal <- list()

    for (gi in seq_along(levels_g)) {
        g <- levels_g[gi]
        cells_g <- names(group_ids)[group_ids == g]
        if (length(cells_g) < min_cells) {
            warning("group '", g, "' skipped: ", length(cells_g),
                    " cells < min_cells = ", min_cells, call. = FALSE)
            next
        }
        if (verbose) message("[OptimizeGammaSuperCell] group '", g, "': ",
                             length(gammas), " granularities")
        counts_g <- .get_counts_matrix(seurat_obj, assay, slot, layer, cells_g)
        obj_g <- subset(seurat_obj, cells = cells_g)
        red <- .sc_prepare_reduction(obj_g, reduction, dims, n_top_genes, n_pc,
                                     seed)

        mem_cols <- list()
        for (k in seq_along(gammas)) {
            gamma <- gammas[k]
            sc <- .sc_run_simplify(obj_g, gamma, red$reduction, red$dims,
                                   assay, k.knn)
            codes <- sc$membership[cells_g]
            uid <- sort(unique(codes))
            lab <- paste0("MC", match(codes, uid))
            names(lab) <- cells_g
            mem_cols[[k]] <- lab
        }
        # gammas collapsing to a single metacell cannot be scored by mcRigor
        n_unique <- vapply(mem_cols, function(l) length(unique(l)), integer(1))
        keep_idx <- which(n_unique >= 2L)
        keep_gamma <- gammas[keep_idx]
        drop_gamma <- gammas[n_unique < 2L]
        if (length(drop_gamma)) {
            warning("group '", g, "': gamma ", paste(drop_gamma, collapse = ", "),
                    " yield(s) a single metacell; mcRigor cannot score it  -  dropped",
                    call. = FALSE)
        }
        if (!length(keep_gamma)) {
            warning("group '", g, "': no candidate granularity produced >= 2 metacells",
                    call. = FALSE)
            next
        }

        code_df <- as.data.frame(mem_cols[keep_idx], stringsAsFactors = FALSE,
                                 check.names = FALSE)
        colnames(code_df) <- as.character(keep_gamma)
        rownames(code_df) <- cells_g

        obj_qc <- Seurat::CreateSeuratObject(
            counts = counts_g, assay = "RNA", min.features = 0, min.cells = 0
        )
        args <- list(Nrep = 50, draw = FALSE, step_save = FALSE,
                     output_file = tempfile(fileext = ".rds"))
        args <- utils::modifyList(args, mcrigor_args)
        res <- tryCatch(
            do.call(
                mcRigor::mcRigor_OPTIMIZE,
                utils::modifyList(
                    list(obj_singlecell = obj_qc, cell_membership = code_df),
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
        best <- .sc_optimize_pick(res, keep_gamma)
        optimal[[g]] <- best$gamma
        sc_vec <- rep(NA_real_, length(gammas))
        sc_vec[match(as.character(best$score_tab$gamma), as.character(gammas))] <-
            best$score_tab$score
        score_rows[[g]] <- data.frame(group = g, gamma = gammas, score = sc_vec)
        # membership labels for the chosen granularity, mapped back to
        # <group>_<i> names in code order
        best_col <- mem_cols[[match(as.character(best$gamma),
                                    as.character(keep_gamma))]]
        best_uid <- sort(unique(best_col))
        g_prefix <- .mc_sanitize(paste0(paste(group.by, collapse = "_"), "_", g))
        membership_cols[[g]] <- paste0(g_prefix, "_", match(best_col, best_uid))
        names(membership_cols[[g]]) <- cells_g
    }

    list(
        optimal_gamma = optimal,
        scores = do.call(rbind, score_rows),
        membership = if (length(membership_cols) > 0) membership_cols else NULL
    )
}

.sc_optimize_pick <- function(res, gammas) {
    # Extract per-gamma scores from mcRigor_OPTIMIZE output (pinned against
    # HEAD c70cd6f: res$scores is a data.frame(gamma, DubRate, ZeroRate,
    # Score); res$best_granularity_level is the winning gamma).
    tab <- if (!is.null(res$scores) && is.data.frame(res$scores)) {
        res$scores
    } else {
        dfs <- Filter(is.data.frame, res)
        if (length(dfs) > 0) dfs[[1]] else NULL
    }
    if (is.null(tab)) {
        stop("unrecognized mcRigor_OPTIMIZE return structure")
    }
    score_col <- intersect(c("Score", "score", "opt_score"), colnames(tab))[1]
    gamma_col <- intersect(c("gamma", "Gamma", "granularity",
                             "best_granularity_level"), colnames(tab))[1]
    if (is.na(score_col) || is.na(gamma_col)) {
        stop("unrecognized mcRigor_OPTIMIZE scores table")
    }
    gvals <- as.numeric(tab[[gamma_col]])
    scores <- as.numeric(tab[[score_col]])
    aligned <- scores[match(as.character(gammas), as.character(gvals))]
    score_tab <- data.frame(gamma = gammas,
                            score = if (!all(is.na(aligned))) aligned else
                              rep(NA_real_, length(gammas)),
                            stringsAsFactors = FALSE)
    best_gamma <- if (!is.null(res$best_granularity_level)) {
        as.numeric(res$best_granularity_level)
    } else if (any(!is.na(score_tab$score))) {
        score_tab$gamma[which.max(score_tab$score)]
    } else {
        NA_real_
    }
    list(gamma = best_gamma, score_tab = score_tab)
}

#' Access metacell results stored by MetacellsByGroupsSuperCell
#'
#' @param seurat_obj A Seurat object processed by
#'   [MetacellsByGroupsSuperCell()].
#' @param wgcna_name Name under which results were stored. Default
#'   `"SuperCell"`.
#'
#' @return The stored list: `metacell_object`, `soft_assignments` (NULL for
#'   this backend), `qc`, `fit_stats`, `zerorate_overall`, `dubrate_overall`,
#'   `score_overall`, `params`, `skipped_groups`, `call`.
#'
#' @export
GetMetacellSuperCell <- function(seurat_obj, wgcna_name = "SuperCell") {
    out <- seurat_obj@misc[[wgcna_name]]
    if (is.null(out)) {
        stop("no metacell results under misc[['", wgcna_name,
             "']]; run MetacellsByGroupsSuperCell first")
    }
    out
}
