# MetacellsByGroupsSEACells: hdWGCNA-style group-stratified metacell
# construction with SEACells as the aggregation backend and mcRigor QC.
#
# Interface mirrors hdWGCNA::MetacellsByGroups (Morabito et al. 2023):
# group.by / ident.group / reduction / dims / assay / slot / layer / mode /
# min_cells / target_metacells / wgcna_name. The aggregation algorithm itself
# is SEACells (Persad et al. 2023) run per group inside a basilisk-managed
# Python environment; quality governance is mcRigor (Liu & Li 2025).

.mc_sanitize <- function(x) {
    x <- gsub("[^[:alnum:]]+", "_", x)
    x <- gsub("_+", "_", x)
    x <- gsub("^_|_$", "", x)
    x
}

.mc_get_group_ids <- function(seurat_obj, group.by) {
    cols <- lapply(group.by, function(col) {
        meta_col <- seurat_obj[[col]][[1]]
        as.character(meta_col)
    })
    ids <- do.call(paste, c(cols, sep = "_"))
    names(ids) <- colnames(seurat_obj)
    ids
}

.mc_n_seacells_for_group <- function(n_cells, cells_per_metacell, target_metacells) {
    n <- floor(n_cells / cells_per_metacell)
    if (!is.null(target_metacells)) {
        n <- min(n, target_metacells)
    }
    max(1, n)
}

.run_seacells_in_basilisk <- function(paths, out_dir, n_seacells, params, seed) {
    helper <- system.file("python", "seacells_helpers.py", package = "cellprogramsMC")
    if (nchar(helper) == 0) {
        stop("seacells_helpers.py not found in installed package")
    }
    proc <- basilisk::basiliskStart(seacells_env, testload = "SEACells")
    on.exit(try(basilisk::basiliskStop(proc), silent = TRUE), add = TRUE)
    out <- basilisk::basiliskRun(
        proc,
        fun = function(helper_path, mtx_path, cells_csv, genes_csv, out_dir,
                       n_seacells, n_top_genes, n_pc, pca_csv,
                       n_waypoint_eigs, waypt_proportion, n_neighbors,
                       convergence_epsilon, max_iter, seed, use_sparse, mode) {
            reticulate::source_python(helper_path, envir = environment())
            run_seacells_group(
                mtx_path = mtx_path, cells_csv = cells_csv, genes_csv = genes_csv,
                out_dir = out_dir, n_seacells = n_seacells,
                n_top_genes = n_top_genes, n_pc = n_pc, pca_csv = pca_csv,
                n_waypoint_eigs = n_waypoint_eigs,
                waypt_proportion = waypt_proportion, n_neighbors = n_neighbors,
                convergence_epsilon = convergence_epsilon, max_iter = max_iter,
                seed = seed, use_sparse = use_sparse, mode = mode
            )
        },
        helper_path = helper,
        mtx_path = paths$mtx, cells_csv = paths$cells, genes_csv = paths$genes,
        out_dir = out_dir, n_seacells = n_seacells,
        n_top_genes = params$n_top_genes, n_pc = params$n_pc,
        pca_csv = paths$pca, n_waypoint_eigs = params$n_waypoint_eigs,
        waypt_proportion = params$waypt_proportion,
        n_neighbors = params$n_neighbors,
        convergence_epsilon = params$convergence_epsilon,
        max_iter = params$max_iter, seed = as.integer(seed),
        use_sparse = params$use_sparse, mode = params$mode
    )
    out
}

#' Construct metacells by groups with SEACells and mcRigor QC
#'
#' Constructs transcriptomically homogeneous metacells from a Seurat object
#' using SEACells (Persad et al. 2023) as the aggregation backend, computed
#' separately within each level of the supplied groupings (typically donor,
#' cell type, or condition), following the hdWGCNA `MetacellsByGroups`
#' interface (Morabito et al. 2023). Metacell quality is governed by mcRigor
#' (Liu & Li 2025): per-metacell mcDiv dubious detection with double-permutation
#' nulls, DubRate/ZeroRate scoring, and optional granularity optimization via
#' [OptimizeGammaSEACells()].
#'
#' The Python SEACells backend runs inside an isolated basilisk environment
#' ([seacells_env]); no Python knowledge is required and the R interface is
#' self-contained. SEACells is GPLv2 and is used strictly as an external
#' dependency, never source-copied.
#'
#' @param seurat_obj A Seurat object with raw counts and a dimensionality
#'   reduction (default `pca`).
#' @param group.by Character vector of metadata column names defining the
#'   stratification (e.g. `c("donor", "celltype")`). Metacells are constructed
#'   separately within each combination of levels.
#' @param ident.group Optional metadata column name used as the metacell-name
#'   prefix; defaults to the interaction of `group.by`.
#' @param reduction Name of the dimensionality reduction to build the SEACells
#'   kernel on. The stored embeddings are exported and used directly
#'   (`build_kernel_on = 'X_pca'`). Set to `NULL` to let the Python side
#'   compute a scanpy PCA (50 PCs on HVGs) instead.
#' @param dims Dimensions of `reduction` to use (indices or names); default
#'   all.
#' @param assay Assay to aggregate. Default `RNA`.
#' @param slot Slot (Seurat v4) to extract counts from. Default `counts`.
#' @param layer Layer (Seurat v5) to extract counts from. Default `counts`.
#' @param mode How to aggregate constituent cells into metacells: `"sum"`
#'   (SEACells default, raw-count sums) or `"average"` (divide by metacell
#'   size). Default `"sum"`.
#' @param cells_per_metacell Target number of single cells per metacell
#'   (SEACells guidance: 75; hdWGCNA k guidance: 20-75). Default 75.
#' @param target_metacells Optional cap on the number of metacells per group.
#' @param min_cells Groups with fewer than this many cells are skipped with a
#'   warning. Default 100.
#' @param n_top_genes Number of highly variable genes used for the scanpy
#'   preprocessing on the Python side. Default 1500.
#' @param n_pc Number of PCs computed by scanpy when `reduction = NULL`.
#'   Default 50.
#' @param n_waypoint_eigs Number of eigenvectors for waypoint initialization.
#'   Default 10.
#' @param waypt_proportion Proportion of archetypes initialized by waypoint
#'   analysis. Default 0.9.
#' @param n_neighbors Number of neighbors for the SEACells graph. Default 15.
#' @param convergence_epsilon Franke-Wolfe convergence threshold. Default 1e-5
#'   (hdWGCNA vignette value).
#' @param max_iter Maximum SEACells fit iterations. Default 100 (the SEACells
#'   0.3.3 fork default; non-convergence is recorded in fit_stats, not fatal).
#' @param use_sparse Use sparse matrix operations in SEACells (CPU only).
#'   Default TRUE.
#' @param use_soft_assignments If TRUE (default), export the cells x metacells
#'   soft-assignment matrix (`A_.T`) and store it on the output.
#' @param run_mcrigor If TRUE (default), run mcRigor dubious-metacell
#'   detection per group and attach DubRate/ZeroRate/Score to the QC table.
#' @param mcrigor_args Optional named list overriding mcRigor_DETECT
#'   arguments: `Nrep` (50), `gene_filter` (0.1), `feature_use` (2000),
#'   `test_cutoff` (0.01), `cor_method` ("pearson"), `aggregate_method`
#'   (defaults to `"sum"`/`"mean"` following `mode`).
#' @param seed Random seed for SEACells initialization. Default 20261016.
#' @param name String appended to metacell names. Default `"SEAC"`.
#' @param wgcna_name Name under which results are stored in
#'   `seurat_obj@misc`. Default `"SEACells"`.
#' @param verbose Print progress. Default TRUE.
#'
#' @return The Seurat object with results stored in `seurat_obj@misc[[wgcna_name]]`:
#'   `metacell_object` (Seurat object of metacells, counts = aggregated
#'   expression, metadata = group, n_cells, mcRigor QC columns),
#'   `soft_assignments` (dgCMatrix, cells x metacells, if
#'   `use_soft_assignments`), `qc` (per-metacell QC data.frame), `params`
#'   (all call parameters), `skipped_groups` (data.frame), and `call`.
#'
#' @details Groups are processed independently; each group's cells are
#'   exported (MatrixMarket + CSV, the hdWGCNA vignette round-trip pattern),
#'   SEACells is fit with `n_SEACells = floor(n_cells / cells_per_metacell)`
#'   (capped by `target_metacells`), metacells are summarized from raw counts,
#'   and mcRigor QC is computed against the group's single-cell counts.
#'   Metacells are named `<group>_<i>` and column-bound across groups. The
#'   result is hdWGCNA-interoperable: pass
#'   `hdWGCNA::SetMetacellObject(obj, GetMetacellSEACells(seurat_obj)$metacell_object)`
#'   after `hdWGCNA::SetupForWGCNA` to use it in a WGCNA network.
#'
#' @references
#' Persad et al. (2023) Nat Biotechnol, doi:10.1038/s41587-023-01716-9;
#' Liu & Li (2025) Nat Commun, doi:10.1038/s41467-025-63626-5;
#' Morabito et al. (2023) Cell Rep Methods, doi:10.1016/j.crmeth.2023.100498.
#'
#' @export
MetacellsByGroupsSEACells <- function(
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
    n_top_genes = 1500,
    n_pc = 50,
    n_waypoint_eigs = 10,
    waypt_proportion = 0.9,
    n_neighbors = 15,
    convergence_epsilon = 1e-5,
    max_iter = 20,
    use_sparse = TRUE,
    use_soft_assignments = TRUE,
    run_mcrigor = TRUE,
    mcrigor_args = list(),
    seed = 20261016,
    name = "SEAC",
    wgcna_name = "SEACells",
    verbose = TRUE
) {
    mode <- match.arg(mode)
    if (!methods::is(seurat_obj, "Seurat")) {
        stop("seurat_obj must be a Seurat object")
    }
    missing_cols <- setdiff(group.by, colnames(seurat_obj@meta.data))
    if (length(missing_cols) > 0) {
        stop("group.by columns not found in meta.data: ",
             paste(missing_cols, collapse = ", "))
    }
    if (!is.null(reduction) && !(reduction %in% names(seurat_obj@reductions))) {
        stop("reduction '", reduction, "' not found; available: ",
             paste(names(seurat_obj@reductions), collapse = ", "))
    }

    group_ids <- .mc_get_group_ids(seurat_obj, group.by)
    levels_g <- sort(unique(group_ids))
    prefix_base <- if (!is.null(ident.group)) ident.group else paste(group.by, collapse = "_")

    params <- list(
        group.by = group.by, ident.group = ident.group, reduction = reduction,
        dims = dims, assay = assay, slot = slot, layer = layer, mode = mode,
        cells_per_metacell = cells_per_metacell,
        target_metacells = target_metacells, min_cells = min_cells,
        n_top_genes = n_top_genes, n_pc = n_pc,
        n_waypoint_eigs = n_waypoint_eigs, waypt_proportion = waypt_proportion,
        n_neighbors = n_neighbors, convergence_epsilon = convergence_epsilon,
        max_iter = max_iter, use_sparse = use_sparse,
        use_soft_assignments = use_soft_assignments,
        run_mcrigor = run_mcrigor, seed = seed, name = name
    )

    mc_list <- list(); obs_list <- list(); soft_list <- list(); qc_list <- list()
    skipped <- list(); fit_stats_list <- list()
    tmp_root <- file.path(tempdir(), paste0("cpmc_", as.integer(Sys.time())))
    dir.create(tmp_root, recursive = TRUE, showWarnings = FALSE)

    for (gi in seq_along(levels_g)) {
        g <- levels_g[gi]
        cells_g <- names(group_ids)[group_ids == g]
        g_prefix <- .mc_sanitize(paste0(prefix_base, "_", g))
        if (verbose) message("[", gi, "/", length(levels_g), "] group '", g,
                             "': ", length(cells_g), " cells")
        if (length(cells_g) < min_cells) {
            skipped[[g]] <- data.frame(group = g, n_cells = length(cells_g),
                                       reason = "below min_cells")
            next
        }

        counts_g <- .get_counts_matrix(seurat_obj, assay, slot, layer, cells_g)
        emb_g <- if (!is.null(reduction)) {
            .get_reduction_embeddings(seurat_obj, reduction, dims, cells_g)
        } else {
            NULL
        }
        out_dir <- file.path(tmp_root, paste0("g", gi))
        paths <- .export_group_for_python(counts_g, emb_g, out_dir)

        n_seacells <- .mc_n_seacells_for_group(
            length(cells_g), cells_per_metacell, target_metacells
        )
        res <- .run_seacells_in_basilisk(paths, out_dir, n_seacells, params, seed)

        mc_out <- .read_metacell_outputs(out_dir)
        mc_mtx <- mc_out$mtx  # genes x metacells
        mc_obs <- mc_out$obs
        m <- ncol(mc_mtx)
        mc_names <- paste0(g_prefix, "_", seq_len(m))
        colnames(mc_mtx) <- mc_names
        rownames(mc_obs) <- mc_names
        mc_obs$Group <- g

        if (!is.null(mc_out$soft) && use_soft_assignments) {
            soft_g <- mc_out$soft  # cells x metacells
            if (nrow(soft_g) == length(cells_g)) {
                rownames(soft_g) <- cells_g
                colnames(soft_g) <- mc_names
                soft_list[[g]] <- soft_g
            } else {
                warning("soft assignment row mismatch in group '", g,
                        "'; dropping soft assignments for this group")
            }
        }

        # ---- mcRigor QC (native R) ----
        qc_g <- NULL
        if (run_mcrigor) {
            qc_g <- tryCatch(
                .mcrigor_detect_group(
                    counts_g = counts_g,
                    membership = .hard_labels_vector(mc_out, cells_g, mc_names),
                    mode = mode, mcrigor_args = mcrigor_args, group_id = g
                ),
                error = function(e) {
                    warning("mcRigor QC failed for group '", g, "': ", conditionMessage(e))
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
        # ZeroRate + Score (mcRigor paper: Score = 1 - w*DubRate - (1-w)*ZeroRate)
        zr_mc <- 1 - Matrix::colSums(mc_mtx != 0) / nrow(mc_mtx)
        qc_g$zerorate_mc <- as.numeric(zr_mc)[match(qc_g$metacell, mc_names)]
        qc_list[[g]] <- qc_g

        mc_list[[g]] <- mc_mtx
        obs_list[[g]] <- mc_obs
        fit_stats_list[[g]] <- data.frame(group = g, mc_out$fit_stats,
                                          stringsAsFactors = FALSE)
    }

    if (length(mc_list) == 0) {
        stop("no groups passed min_cells = ", min_cells)
    }

    # cbind2 is a BINARY S4 generic: do.call with >2 groups passes N args and
    # the selected method fails ('argument "x" is missing'). Reduce instead.
    mc_all <- Reduce(Matrix::cbind2, mc_list)
    obs_all <- do.call(rbind, obs_list)
    qc_all <- do.call(rbind, qc_list)
    soft_all <- if (length(soft_list) > 0 && use_soft_assignments) {
        do.call(rbind, soft_list)  # cells x metacells
    } else {
        NULL
    }
    fit_stats_all <- do.call(rbind, fit_stats_list)

    # overall ZeroRate of the metacell matrix
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

    skipped_df <- if (length(skipped) > 0) do.call(rbind, skipped) else
        data.frame(group = character(0), n_cells = integer(0), reason = character(0))

    out <- list(
        metacell_object = m_obj,
        soft_assignments = soft_all,
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

.hard_labels_vector <- function(mc_out, cells_g, mc_names) {
    # Reconstruct the per-cell hard metacell label vector from the soft
    # assignment matrix (row-wise argmax over the sparse matrix), used as
    # mcRigor's cell_membership input.
    soft <- mc_out$soft
    if (is.null(soft)) {
        stop("soft assignments required for mcRigor membership reconstruction")
    }
    # the soft matrix is exported as a bare MTX (no dimnames); rows are in the
    # exported cell order (ad.obs order = cells_g order)
    if (is.null(rownames(soft))) {
        if (nrow(soft) != length(cells_g)) {
            stop("soft assignment rows do not match group cells")
        }
        rownames(soft) <- cells_g
    }
    s <- as(soft, "dgTMatrix")
    trip <- data.frame(i = s@i + 1L, j = s@j + 1L, x = s@x)
    # row-wise argmax: order by row, then descending weight
    ord <- order(trip$i, -trip$x)
    trip <- trip[ord, ]
    first <- !duplicated(trip$i)
    lab_idx <- trip$j[first]
    lab_row <- trip$i[first]
    lab <- rep(NA_character_, nrow(soft))
    lab[lab_row] <- mc_names[lab_idx]
    names(lab) <- rownames(soft)
    lab[is.na(lab)] <- mc_names[1]  # all-zero rows: argmax undefined -> first
    lab[cells_g]
}
