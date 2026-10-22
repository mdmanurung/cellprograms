# Covariate prediction is a trained transformation, not model.matrix(query).
# Only these documented deterministic transforms are accepted. Custom functions
# may hide query-dependent state and must be applied outside the fit explicitly.
.validate_covariate_formula <- function(covariates) {
  f <- if (inherits(covariates, "formula")) covariates else stats::reformulate(covariates)
  if (length(f) != 2L) cli::cli_abort("Covariates require a one-sided formula.")
  walk <- function(e) {
    if (!is.call(e)) return(invisible(NULL))
    head <- e[[1L]]
    name <- if (is.call(head) && identical(head[[1L]], as.name("::"))) {
      paste0(head[[2L]], "::", head[[3L]])
    } else as.character(head)
    allowed <- c("~", "+", "-", "*", "/", "^", ":", "(", "I", "log", "log1p",
                 "sqrt", "poly", "stats::poly", "ns", "splines::ns", "bs",
                 "splines::bs", "scale")
    if (length(name) != 1L || !name %in% allowed) {
      cli::cli_abort("Unsupported covariate transform: {.val {paste(name, collapse='')}}. Precompute custom transforms with an explicit trained recipe before fitting.")
    }
    if (name == "scale" && (length(e) != 2L || !is.symbol(e[[2L]]))) {
      cli::cli_abort("Only scale(variable) is supported in covariate formulas; precompute other scaling recipes.")
    }
    for (a in as.list(e)[-1L]) walk(a)
    invisible(NULL)
  }
  walk(f)
  invisible(f)
}

.trained_covariate_formula <- function(f, m) {
  .validate_covariate_formula(f)
  walk <- function(e) {
    if (!is.call(e)) return(e)
    if (identical(e[[1L]], as.name("scale"))) {
      z <- m[[as.character(e[[2L]])]]
      if (!is.numeric(z) || !is.finite(stats::sd(z)) || stats::sd(z) <= 0) {
        cli::cli_abort("scale(variable) needs a varying numeric training covariate.")
      }
      return(substitute(scale(V, center = C, scale = S),
                        list(V = e[[2L]], C = mean(z), S = stats::sd(z))))
    }
    as.call(lapply(as.list(e), walk))
  }
  f[[2L]] <- walk(f[[2L]])
  # Do not serialize a caller's mutable workspace as prediction state.
  environment(f) <- list2env(list(poly = stats::poly, ns = splines::ns,
                                  bs = splines::bs), parent = baseenv())
  f
}

.predict_covariate_design <- function(state, donors, meta, ct) {
  if (!is.list(state) || !identical(state$schema_version, "1.0.0")) {
    cli::cli_abort("Missing trained covariate design state for {.val {ct}}; refit legacy residualized models before projection.")
  }
  if (anyDuplicated(meta$observation_id) || anyNA(meta$observation_id)) {
    cli::cli_abort("Query sample metadata must have unique nonmissing observation_id values.")
  }
  absent <- setdiff(unique(donors), meta$observation_id)
  if (length(absent)) cli::cli_abort("Query donor(s) missing from sample_metadata: {.val {absent}}.")
  if (length(setdiff(state$variables, names(meta)))) {
    cli::cli_abort("Query sample_metadata lacks trained covariates.")
  }
  m <- meta[match(donors, meta$observation_id), state$variables, drop = FALSE]
  if (anyNA(m)) cli::cli_abort("Missing covariate values for query donors of cell type {.val {ct}}.")
  # model.frame uses the training predvars (poly coefficients / spline knots)
  # and training factor levels, even for a single/repeated query row.
  mf <- stats::model.frame(state$terms, m, xlev = state$xlevels,
                           na.action = stats::na.fail)
  X <- stats::model.matrix(state$terms, mf, contrasts.arg = state$contrasts)
  if (!identical(colnames(X), state$columns) || any(!is.finite(X))) {
    cli::cli_abort("Replayed trained covariate design does not match the recorded finite columns.")
  }
  X
}
