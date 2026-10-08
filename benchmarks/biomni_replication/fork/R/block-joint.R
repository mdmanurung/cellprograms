## block-joint.R — Level 3 candidate A: block-structured joint EBMF
##
## Greedy two-pass strategy built from flashier primitives:
##   shared_first:  fit EBMF on the stacked matrix -> classify factors by
##                  block-loading concentration -> backfit each cell type on
##                  its residual to recover private structure.
##   private_first: take residuals of the stage-1 per-cell-type fits, stack,
##                  and fit shared structure on the stacked residuals.
##
## Levels:
##   expression: stacked matrix rows = union of observations, columns =
##               cell-type-prefixed genes (NA for missing sample x cell type).
##   scores:     the combined program-score matrix Z_all (column blocks =
##               cell types), each block scaled to unit Frobenius norm.

.build_stacked <- function(fit, level = c("expression", "scores")) {
  level <- match.arg(level)
  cts <- names(fit$input$matrices)
  ids <- unique(unlist(fit$input$observation_ids))

  if (level == "expression") {
    blocks <- fit$input$matrices
  } else {
    blocks <- fit$scores
    blocks <- blocks[vapply(blocks, ncol, integer(1L)) > 0L]
    cts <- names(blocks)
  }

  ncols <- vapply(blocks, ncol, integer(1L))
  X <- matrix(NA_real_, nrow = length(ids), ncol = sum(ncols),
              dimnames = list(ids, unlist(lapply(cts, function(ct) {
                paste0(ct, "__", colnames(blocks[[ct]]))
              }))))
  block_id <- rep(cts, times = ncols)
  for (ct in cts) {
    B <- blocks[[ct]]
    X[rownames(B), block_id == ct] <- B
  }
  block_scales <- rep(1, length(cts)); names(block_scales) <- cts
  if (level == "scores") {
    ## unit-Frobenius block scaling so no cell type dominates by scale
    for (ct in cts) {
      idx <- block_id == ct
      frob <- sqrt(sum(X[, idx]^2, na.rm = TRUE))
      if (frob > 0) {
        X[, idx] <- X[, idx] / frob
        block_scales[ct] <- frob
      }
    }
  }
  list(X = X, block_id = block_id, block_scales = block_scales, level = level)
}

.classify_factors <- function(W, block_id, private_threshold = 0.8,
                              global_entropy = 0.8) {
  K <- ncol(W)
  cts <- unique(block_id)
  C <- length(cts)
  out <- vector("list", K)
  for (k in seq_len(K)) {
    w2 <- W[, k]^2
    mass <- vapply(cts, function(ct) sum(w2[block_id == ct]), numeric(1L))
    mass <- mass / max(sum(mass), 1e-12)
    H <- -sum(ifelse(mass > 0, mass * log(mass), 0)) / log(C)
    cls <- if (max(mass) > private_threshold) {
      paste0("private:", cts[which.max(mass)])
    } else if (C > 2 && H > global_entropy) {
      "global"
    } else {
      paste0("partial:", paste(cts[mass > 0.1], collapse = "+"))
    }
    out[[k]] <- list(class = cls, block_mass = mass, entropy = H)
  }
  out
}

#' Block-structured joint EBMF (Level 3 candidate A)
#'
#' @param fit A `cell_program_fit` (uses its processed matrices / scores).
#' @param level `"expression"` or `"scores"`.
#' @param order `"shared_first"` (default) or `"private_first"` (ablation).
#' @param max_factors Candidate factors for the joint pass.
#' @param max_private Candidate factors per cell type in the private pass.
#' @param private_threshold Block-mass threshold for private classification.
#' @param global_entropy Normalized-entropy threshold for global classification.
#' @param seed Random seed.
#' @param verbose Verbosity.
#' @return An integration-style list with method `"block"`.
#' @export
fit_joint_block <- function(fit, level = c("scores", "expression"),
                            order = c("shared_first", "private_first"),
                            max_factors = 20L, max_private = 10L,
                            private_threshold = 0.8, global_entropy = 0.8,
                            seed = NULL, verbose = 0L) {
  level <- match.arg(level)
  order <- match.arg(order)
  if (!is.null(seed)) set.seed(seed)
  st <- .build_stacked(fit, level)
  X <- st$X; block_id <- st$block_id
  ebnm_list <- list(ebnm::ebnm_normal, ebnm::ebnm_point_normal)

  shared_fit <- NULL; private_fits <- NULL

  if (order == "shared_first") {
    shared_fit <- flashier::flash(X, ebnm_fn = ebnm_list, var_type = 2L,
                                  greedy_Kmax = max_factors, backfit = TRUE,
                                  nullcheck = TRUE, verbose = verbose)
    ## private pass: per-cell-type residuals after removing shared reconstruction.
    ## Expression level only: at scores level each block has too few columns
    ## for EB prior estimation, and private structure IS the stage-1 programs.
    private_fits <- list()
    if (shared_fit$n_factors > 0L) {
      recon <- shared_fit$L_pm %*% t(shared_fit$F_pm)
    } else {
      recon <- matrix(0, nrow(X), ncol(X))
    }
    if (level == "scores") {
      private_fits <- NULL
    }
    for (ct in if (level == "expression") unique(block_id) else character(0)) {
      idx <- which(block_id == ct)
      obs <- rownames(fit$input$matrices[[ct]])
      if (level == "expression") {
        Yc <- X[obs, idx, drop = FALSE] - recon[obs, idx, drop = FALSE]
      } else {
        Yc <- X[, idx, drop = FALSE] - recon[, idx, drop = FALSE]
        Yc <- Yc[stats::complete.cases(Yc), , drop = FALSE]
      }
      if (nrow(Yc) < 5L) next
      ## skip blocks whose residual is numerically zero (fully explained by
      ## the shared pass): EBNM standard errors degenerate on constant input
      resid_rms <- sqrt(mean(Yc^2, na.rm = TRUE))
      if (!is.finite(resid_rms) || resid_rms < 1e-8) next
      private_fits[[ct]] <- tryCatch(
        flashier::flash(Yc, ebnm_fn = ebnm_list, var_type = 2L,
                        greedy_Kmax = max_private, backfit = TRUE,
                        nullcheck = TRUE, verbose = verbose),
        error = function(e) {
          .warnf("block private pass failed for '%s': %s", ct, conditionMessage(e))
          NULL
        })
    }
  } else {
    ## private_first: residuals of stage-1 fits, stacked, then shared pass
    if (level == "expression") {
      R <- X
      for (ct in names(fit$fits)) {
        idx <- which(block_id == ct)
        obs <- rownames(fit$input$matrices[[ct]])
        K <- fit$fits[[ct]]$n_factors
        if (K > 0L) {
          R[obs, idx] <- fit$input$matrices[[ct]] -
            fit$scores[[ct]] %*% t(fit$loadings[[ct]])
        }
      }
    } else {
      R <- X  ## stage-1 scores are already "private" factors; shared pass on residuals
      ## residualize each block against its own column mean structure is a no-op;
      ## the shared pass operates on the raw score blocks
    }
    shared_fit <- flashier::flash(R, ebnm_fn = ebnm_list, var_type = 2L,
                                  greedy_Kmax = max_factors, backfit = TRUE,
                                  nullcheck = TRUE, verbose = verbose)
  }

  ## classify shared factors
  cls <- list()
  if (!is.null(shared_fit) && shared_fit$n_factors > 0L) {
    cls <- .classify_factors(shared_fit$F_pm, block_id,
                             private_threshold, global_entropy)
  }

  summ <- data.frame(
    factor_id = sprintf("SF%02d", seq_len(shared_fit$n_factors %||% 0L)),
    class = vapply(cls, `[[`, character(1L), "class"),
    entropy = vapply(cls, `[[`, numeric(1L), "entropy"),
    stringsAsFactors = FALSE
  )
  if (nrow(summ) > 0L) {
    mass <- do.call(rbind, lapply(cls, `[[`, "block_mass"))
    colnames(mass) <- paste0("mass_", unique(block_id))
    summ <- cbind(summ, mass)
  }

  out <- list(
    method = "block",
    level = level,
    order = order,
    shared_fit = shared_fit,
    private_fits = private_fits,
    block_id = block_id,
    block_scales = st$block_scales,
    summary = summ
  )
  ## expose program weights for trace_state when operating on scores
  if (level == "scores" && shared_fit$n_factors > 0L) {
    pw <- shared_fit$F_pm
    rownames(pw) <- colnames(X)
    colnames(pw) <- summ$factor_id
    out$program_weights <- pw
    out$states <- shared_fit$L_pm
    colnames(out$states) <- summ$factor_id
  }
  if (level == "expression" && shared_fit$n_factors > 0L) {
    gl <- list()
    for (ct in unique(block_id)) {
      W <- shared_fit$F_pm[block_id == ct, , drop = FALSE]
      rownames(W) <- sub(paste0("^", ct, "__"), "", rownames(W))
      colnames(W) <- summ$factor_id
      gl[[ct]] <- W
    }
    out$gene_loadings <- gl
    out$states <- shared_fit$L_pm
    colnames(out$states) <- summ$factor_id
  }
  out
}
