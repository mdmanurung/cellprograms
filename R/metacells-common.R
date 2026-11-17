# Shared helpers for the SuperCell 2.0 + mcRigor metacell path.
#
# mcRigor QC: Liu & Li (2025) Nat Commun, doi:10.1038/s41467-025-63626-5.
# - mcDiv: divergence of the within-metacell feature-correlation matrix from
#   identity; null from a double permutation (within-cell preserving library
#   sizes, then within-feature); size-specific thresholds.
# - DubRate: fraction of dubious metacells; ZeroRate: fraction of zeros in the
#   metacell matrix; Score = 1 - w*DubRate - (1-w)*ZeroRate (w = 0.5 default).


.mc_sanitize <- function(x) {
    x <- gsub("[^[:alnum:]]+", "_", x)
    x <- gsub("_+", "_", x)
    x <- gsub("^_|_$", "", x)
    x
}

.mc_get_group_ids <- function(seurat_obj, group.by) {
    cols <- lapply(group.by, function(col) {
        # Seurat's [[ errors on unknown columns; normalize to a clean abort.
        meta_col <- tryCatch(seurat_obj[[col]][[1]], error = function(e) NULL)
        if (is.null(meta_col)) {
            cli::cli_abort("group.by column {.val {col}} not found in the Seurat object's metadata.")
        }
        as.character(meta_col)
    })
    # paste() renders NA as the literal string "NA", silently lumping
    # unassigned cells into a group named "NA" - fail closed instead.
    na_cells <- Reduce(`|`, lapply(cols, is.na))
    if (any(na_cells)) {
        bad <- group.by[vapply(cols, anyNA, logical(1))]
        cli::cli_abort(
            "{sum(na_cells)} cell(s) have NA in group.by column{?s} {.val {bad}}; subset or fill them before metacell construction."
        )
    }
    ids <- do.call(paste, c(cols, sep = "_"))
    names(ids) <- colnames(seurat_obj)
    ids
}

.mcrigor_detect_group <- function(counts_g, membership, mode, mcrigor_args,
                                  group_id) {
    if (!requireNamespace("mcRigor", quietly = TRUE)) {
        cli::cli_abort(
            "{.pkg mcRigor} is required for metacell QC. Install it with {.code remotes::install_github('JSB-UCLA/mcRigor')}."
        )
    }
    # mcRigor calling convention, pinned against mcRigor HEAD c70cd6f (the
    # benchmark's verbatim-verified construction; the OLDER integer-codes
    # convention broke when DETECT's `named_membership != ""` filter began
    # coercing integer labels to NA -> empty assigned_cellid -> "Cannot find
    # cells provided"). Labels must be STRINGS that do NOT start with a digit:
    # digit-leading strings trigger Seurat's AggregateExpression "g"-prefix
    # mangling and crash buildmc ("No cell overlap between new meta data and
    # Seurat object"). "MC<n>" strings pass both paths. tgamma is left at
    # DEFAULT (passing tgamma="1" re-triggers the integer-coercion crash).
    # NOTE: as.character() drops vector names in R 4.4  -  build the labels
    # directly and reassign names explicitly.
    lab_in <- as.character(membership)
    lab <- paste0("MC", as.integer(factor(lab_in, levels = sort(unique(lab_in)))))
    names(lab) <- names(membership)
    # minimal Seurat object holding the group's single-cell counts
    obj_g <- Seurat::CreateSeuratObject(
        counts = counts_g, assay = "RNA", min.features = 0, min.cells = 0
    )
    membership_df <- data.frame(
        X = lab, stringsAsFactors = FALSE, check.names = FALSE,
        row.names = names(lab)
    )
    colnames(membership_df) <- "1"
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
        cell_membership = membership_df
    )
    res <- do.call(mcRigor::mcRigor_DETECT, utils::modifyList(base_args, args))
    out <- .mcrigor_normalize(res, group_id = group_id)
    # remap MC codes back to the original metacell labels: TabMC rownames are
    # the "MC<n>" labels; code n corresponds to the n-th sorted unique input
    # label (the factor level order used to build the codes)
    code_to_lab <- sort(unique(lab_in))
    code_idx <- suppressWarnings(as.integer(sub("^MC", "", out$metacell)))
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

.mc_validate_gammas <- function(gammas) {
    # OptimizeGammaSuperCell candidate granularities: positive integers,
    # deduplicated and sorted. Duplicates would produce identically named
    # membership columns for mcRigor_OPTIMIZE; non-integers fail cryptically
    # inside SCimplify.
    if (!is.numeric(gammas) || !length(gammas) || any(!is.finite(gammas)) ||
        any(gammas < 1) || any(gammas != round(gammas))) {
        cli::cli_abort("{.arg gammas} must be positive integers, got {.val {gammas}}.")
    }
    sort(unique(as.integer(gammas)))
}

.get_counts_matrix <- function(seurat_obj, assay, slot, layer, cells) {
    # Seurat v5 uses layers; v4 uses slots. Support both.
    counts <- tryCatch(
        {
            lays <- SeuratObject::Layers(seurat_obj[[assay]])
            if (length(lays) > 0) {
                lay <- layer
                if (!(lay %in% lays)) {
                    lay <- grep("counts", lays, value = TRUE)[1]
                    if (is.na(lay)) {
                        cli::cli_abort("Assay {.val {assay}} has no layer {.val {layer}} and no counts-like fallback layer; available layers: {.val {lays}}.")
                    }
                    cli::cli_warn("Layer {.val {layer}} not found in assay {.val {assay}}; falling back to layer {.val {lay}}.")
                }
                SeuratObject::LayerData(seurat_obj[[assay]], layer = lay)
            } else {
                NULL
            }
        },
        error = function(e) NULL
    )
    if (is.null(counts)) {
        counts <- SeuratObject::GetAssayData(seurat_obj[[assay]], slot = slot)
    }
    counts[, cells, drop = FALSE]
}
