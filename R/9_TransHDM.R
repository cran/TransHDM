#' TransHDM: High-Dimensional Mediation Analysis with Transfer Learning
#'
#' @description
#' The \code{TransHDM} function performs high-dimensional mediation analysis
#' under a transfer learning framework. It identifies and estimates indirect
#' (mediation) effects of a high-dimensional set of mediators between an
#' exposure and an outcome by integrating a target dataset with one or more
#' source datasets. Multiple source datasets are row-bound into a single
#' source before screening and estimation.
#'
#' @details
#' \code{TransHDM} implements a three-step workflow that can also be executed
#' step-by-step using the individual components:
#' \enumerate{
#'   \item \strong{Step 1 — Screening:} \code{\link{SIS}} performs Sure
#'     Independence Screening to reduce dimensionality by ranking mediators
#'     based on their marginal associations.
#'   \item \strong{Step 2 — De-biased estimation:}
#'     \code{\link{mediation_inference}} fits de-biased lasso models to
#'     estimate exposure-mediator (alpha) and mediator-outcome (beta) effects.
#'   \item \strong{Step 3 — Multiple testing:} \code{\link{joint_test}}
#'     applies joint multiple testing with FDR control to identify
#'     significant mediators and computes effect summaries.
#' }
#' Set \code{use_SIS = FALSE} to skip Step 1 and use all mediators directly.
#' A list of source datasets is row-bound into one source before these steps.
#' Estimation then proceeds exactly as in the single-source procedure.
#'
#' The variance estimator follows Algorithm 3 of Tian and Feng (2023).
#' For Gaussian linear models, it assumes unit error variance.
#' Under the required source-target heterogeneity conditions, pooled
#' source and target design information is used to estimate the
#' correlation structure, while the debiasing correction is based
#' on target-domain data. This assumption is not automatically
#' verified by the function.
#'
#' @seealso \code{\link{SIS}}, \code{\link{mediation_inference}},
#'   \code{\link{joint_test}}
#'
#' @param target_data A data frame containing the target dataset. All variables must be numeric.
#' @param source_data A data frame, or a list of data frames, containing source
#'   datasets (optional, default: NULL). All variables must be numeric and use
#'   the same column names as \code{target_data}. A list is row-bound into one
#'   source dataset before screening and estimation.
#' @param Y Character string specifying the outcome variable name.
#' @param D Character string specifying the exposure (treatment) variable name.
#' @param M Character vector specifying mediator variable names.
#' @param X Character vector specifying covariate variable names..
#' @param transfer A logical value (default: \code{FALSE}) indicating whether
#'   to enable transfer learning by incorporating information from
#'   \code{source_data}.
#' @param verbose A logical value (default: \code{TRUE}) controlling whether
#'   progress messages are printed to the console.
#' @param ncore An integer (default: 1) specifying the number of CPU cores to
#'   use for parallel computation when fitting mediator models.
#' @param topN An integer (default: \code{NULL}) specifying the number of
#'   mediators to retain after Sure Independence Screening (SIS).
#'   If \code{NULL}, the number is determined automatically based on
#'   the data dimensions.
#' @param use_SIS A logical value (default: \code{TRUE}) indicating whether to apply
#'   Sure Independence Screening (SIS) for initial variable screening. When \code{FALSE},
#'    the full set of mediators is used in estimation and inference.
#' @param p_cutoff A numeric value (default: 0.05) specifying the maximum p-value threshold for
#'   determining significant mediators.
#'
#' @return A list with the following components:
#' \itemize{
#'   \item \code{contributions}: A data frame of identified mediators containing:
#'     \itemize{
#'       \item \code{mediator}: Mediator name
#'       \item \code{alpha}: Estimated exposure mediator effect
#'       \item \code{alpha_pv}: P-value for the exposure mediator effect
#'       \item \code{beta}: Estimated mediator outcome effect
#'       \item \code{beta_pv}: P-value for the mediator outcome effect
#'       \item \code{alpha_beta}: Estimated indirect (mediation) effect
#'       \item \code{ab_pv}: Joint p-value, \code{max(alpha_pv, beta_pv)}
#'       \item \code{fdr_est}: Pointwise FDR from \code{HDMT::fdr_est}
#'       \item \code{pa}: Proportion of the total effect mediated
#'     }
#'   \item \code{effects}: A data frame summarizing the total indirect effect,
#'     direct effect, total effect, and proportion mediated.
#'   \item \code{IDE_est}: A numeric vector of indirect effect estimates for all
#'     specified mediators, with non-selected mediators set to zero.
#'   \item \code{DE_est}: The estimated direct effect of the exposure on the
#'     outcome.
#'   \item \code{all_mediators}: A data frame of all SIS-screened mediators
#'     (including non-significant ones) with columns \code{mediator},
#'     \code{alpha}, \code{alpha_pv}, \code{beta}, \code{beta_pv}, and
#'     \code{significant} (logical).
#' }
#'
#' @references
#' Pan L, Liu Y, Huang C, Lin R, Yu Y, Qin G.
#' Transfer learning reveals the mediating mechanisms of cross-ethnic lipid metabolic
#' pathways in the association between APOE gene and Alzheimer's disease.
#' Brief Bioinform. 2025;26(5):bbaf460.
#' \doi{10.1093/bib/bbaf460}
#'
#' @examples
#' \donttest{
#' set.seed(123)
#'
#' # Target data
#' target_data <- gen_simData_homo(n = 50, p_x = 3, p_m = 20, rho = 0.1)$data
#'
#' # Source data
#' source_data <- gen_simData_homo(n = 100, p_x = 3, p_m = 20, rho = 0.1, source = TRUE,
#' transferable = TRUE)$data
#'
#' # Run TransHDM
#' result <- TransHDM(
#'   target_data = target_data,
#'   source_data = source_data,
#'   Y = "Y",
#'   D = "D",
#'   M = paste0("M", 1:20),
#'   X = paste0("X", 1:3),
#'   transfer = TRUE,
#'   ncore = 1,
#'   topN = 10
#' )
#' summary(result)
#' }
#'
#' @export
#'

TransHDM <- function(
    target_data,
    source_data = NULL,
    Y,
    D,
    M,
    X,
    transfer = FALSE,
    verbose = TRUE,
    ncore = 1,
    topN = NULL,
    use_SIS =TRUE,
    p_cutoff = 0.05
) {

  checked <- .check_analysis_input(
    target_data = target_data,
    source_data = source_data,
    Y = Y, D = D, M = M, X = X,
    transfer = transfer,
    ncore = ncore,
    topN = topN,
    p_cutoff = p_cutoff,
    verbose = verbose,
    use_SIS = use_SIS
  )
  Y <- checked$Y
  D <- checked$D
  M <- checked$M
  X <- checked$X

  source_data <- .combine_sources(
    source_data,
    verbose = verbose && isTRUE(transfer)
  )

  if (use_SIS) {
    screen_result <- SIS(
      target_data = target_data,
      source_data = source_data,
      Y = Y,
      D = D,
      M = M,
      X = X,
      topN = topN,
      transfer = transfer,
      verbose = verbose,
      ncore = ncore
    )
    infer <- mediation_inference(
      screen_result = screen_result,
      transfer = transfer,
      verbose = verbose,
      ncore = ncore
    )
  } else {
    if (verbose) {
      message("Step 1: SIS skipped; using full set of mediators.")
    }
    infer <- mediation_inference(
      target_data = target_data,
      source_data = source_data,
      Y = Y,
      D = D,
      M = M,
      X = X,
      transfer = transfer,
      verbose = verbose,
      ncore = ncore
    )
  }

  joint_test(
    inference_result = infer,
    p_cutoff = p_cutoff,
    verbose = verbose
  )
}



.start_cluster <- function(ncore) {
  cl <- parallel::makeCluster(ncore)
  lib_paths <- .libPaths()
  tryCatch(
    {
      parallel::clusterExport(cl, "lib_paths", envir = environment())
      parallel::clusterEvalQ(cl, {
        .libPaths(lib_paths)
        loadNamespace("TransHDM")
        NULL
      })
    },
    error = function(e) {
      parallel::stopCluster(cl)
      stop(e)
    }
  )
  cl
}

.stop_cluster <- function(cl) {
  if (is.null(cl)) {
    return(invisible(NULL))
  }
  parallel::stopCluster(cl)
  invisible(NULL)
}

.combine_sources <- function(source_data, verbose = FALSE) {
  if (is.null(source_data)) {
    return(NULL)
  }
  if (is.data.frame(source_data)) {
    return(source_data)
  }
  if (!is.list(source_data) || length(source_data) == 0) {
    stop("'source_data' must be a data frame or a list of data frames.")
  }
  if (!all(vapply(source_data, is.data.frame, logical(1)))) {
    stop("Each element of 'source_data' must be a data frame.")
  }
  if (length(source_data) == 1) {
    return(source_data[[1]])
  }

  ref <- colnames(source_data[[1]])
  same_cols <- vapply(source_data, function(df) {
    identical(sort(colnames(df)), sort(ref))
  }, logical(1))
  if (!all(same_cols)) {
    stop(
      "All source datasets must contain the same columns. ",
      "Mismatched source index: ",
      paste(which(!same_cols), collapse = ", "),
      "."
    )
  }

  source_data <- lapply(source_data, function(df) df[, ref, drop = FALSE])
  combined <- do.call(rbind, c(source_data, list(make.row.names = FALSE)))

  if (verbose) {
    message(
      "Combining ", length(source_data),
      " source datasets into one source (n = ", nrow(combined), ")."
    )
  }
  combined
}

.resolve_cols <- function(data, cols, name) {
  if (is.null(cols)) return(NULL)

  if (is.numeric(cols)) {
    if (any(cols < 1 | cols > ncol(data)))
      stop(sprintf("%s index out of range", name))
    return(colnames(data)[cols])
  }

  if (is.character(cols)) {
    if (!all(cols %in% colnames(data))) {
      stop(sprintf(
        "Some %s variables not found in data: %s",
        name,
        paste(setdiff(cols, colnames(data)), collapse = ", ")
      ))
    }
    return(cols)
  }

  stop(sprintf("%s must be column names or column indices", name))
}

.check_analysis_input <- function(
    target_data,
    source_data = NULL,
    Y, D, M, X = NULL,
    transfer = FALSE,
    ncore = 1,
    topN = NULL,
    p_cutoff = NULL,
    kfold = NULL,
    t = NULL,
    verbose = TRUE,
    use_SIS = NULL,
    require_source = FALSE
) {
  if (!is.data.frame(target_data) || nrow(target_data) < 2) {
    stop("'target_data' must be a data frame with at least 2 rows.")
  }

  .check_flag <- function(val, name) {
    if (is.null(val)) return(invisible(NULL))
    if (!(is.logical(val) && length(val) == 1L && !is.na(val))) {
      stop(sprintf("'%s' must be TRUE or FALSE.", name))
    }
  }
  .check_flag(transfer, "transfer")
  .check_flag(verbose, "verbose")
  .check_flag(use_SIS, "use_SIS")

  if (!is.numeric(ncore) || length(ncore) != 1L || is.na(ncore) ||
      ncore < 1 || ncore != as.integer(ncore)) {
    stop("'ncore' must be a positive integer.")
  }
  if (!is.null(topN) &&
      (!is.numeric(topN) || length(topN) != 1L || is.na(topN) ||
       topN < 1 || topN != as.integer(topN))) {
    stop("'topN' must be NULL or a positive integer.")
  }
  if (!is.null(p_cutoff) &&
      (!is.numeric(p_cutoff) || length(p_cutoff) != 1L || is.na(p_cutoff) ||
       p_cutoff <= 0 || p_cutoff > 1)) {
    stop("'p_cutoff' must be a single number in (0, 1].")
  }
  if (!is.null(kfold) &&
      (!is.numeric(kfold) || length(kfold) != 1L || is.na(kfold) ||
       kfold < 2 || kfold != as.integer(kfold) || kfold > nrow(target_data))) {
    stop("'kfold' must be an integer >= 2 and no larger than the number of target rows.")
  }
  if (!is.null(t) &&
      (!is.numeric(t) || length(t) != 1L || is.na(t) || t < 0)) {
    stop("'t' must be a single nonnegative number.")
  }

  Y <- .resolve_cols(target_data, Y, "Y")
  D <- .resolve_cols(target_data, D, "D")
  M <- .resolve_cols(target_data, M, "M")
  X <- .resolve_cols(target_data, X, "X")
  if (length(Y) != 1L) stop("'Y' must name exactly one outcome column.")
  if (length(D) != 1L) stop("'D' must name exactly one exposure column.")
  if (length(M) < 1L) stop("At least one mediator must be supplied in 'M'.")

  vars <- c(Y, D, M, X)
  if (anyDuplicated(vars)) {
    stop(
      "'Y', 'D', 'M', and 'X' must name distinct columns. Duplicates: ",
      paste(unique(vars[duplicated(vars)]), collapse = ", ")
    )
  }

  .check_block <- function(data, label) {
    block <- data[, vars, drop = FALSE]
    nonnumeric <- vars[!vapply(block, is.numeric, logical(1))]
    if (length(nonnumeric)) {
      stop(sprintf(
        "%s contains nonnumeric columns: %s",
        label, paste(nonnumeric, collapse = ", ")
      ))
    }
    missing <- vars[vapply(block, anyNA, logical(1))]
    if (length(missing)) {
      stop(sprintf(
        "%s contains missing values in: %s",
        label, paste(missing, collapse = ", ")
      ))
    }
    nonfinite <- vars[vapply(block, function(z) any(is.infinite(z)), logical(1))]
    if (length(nonfinite)) {
      stop(sprintf(
        "%s contains nonfinite values in: %s",
        label, paste(nonfinite, collapse = ", ")
      ))
    }
    predictors <- c(D, M, X)
    zero_var <- predictors[vapply(data[, predictors, drop = FALSE], function(z) {
      stats::var(z) == 0
    }, logical(1))]
    if (length(zero_var)) {
      stop(sprintf(
        "%s contains zero-variance predictors: %s",
        label, paste(zero_var, collapse = ", ")
      ))
    }
    if (stats::var(data[[Y]]) == 0) {
      stop(sprintf("%s outcome '%s' has zero variance.", label, Y))
    }
  }

  .check_block(target_data, "Target data")

  sources <- if (is.null(source_data)) {
    NULL
  } else if (is.data.frame(source_data)) {
    list(source_data)
  } else if (is.list(source_data) && length(source_data) > 0 &&
             all(vapply(source_data, is.data.frame, logical(1)))) {
    source_data
  } else {
    stop("'source_data' must be a data frame or a non-empty list of data frames.")
  }

  if (isTRUE(require_source) && (is.null(sources) || length(sources) == 0)) {
    stop("'source_data' must be a data frame or a non-empty list of source data frames.")
  }
  if (isTRUE(transfer) && (is.null(sources) || length(sources) == 0)) {
    stop("'transfer = TRUE' requires at least one source dataset in 'source_data'.")
  }

  if (!is.null(sources)) {
    needed <- vars
    for (i in seq_along(sources)) {
      missing_cols <- setdiff(needed, colnames(sources[[i]]))
      if (length(missing_cols)) {
        stop(sprintf(
          "Source dataset %d does not contain the same analysis variables as the target data. Missing: %s",
          i, paste(missing_cols, collapse = ", ")
        ))
      }
      .check_block(sources[[i]], sprintf("Source dataset %d", i))
    }
  }

  list(Y = Y, D = D, M = M, X = X)
}

