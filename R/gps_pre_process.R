#' Generalized Propensity Score (GPS) Model Fitting
#'
#' This function estimates generalized propensity scores (GPS) for multi-level
#' or continuous treatment settings using various machine learning or regression
#' models. Supported GPS models include multinomial logistic regression (`logit`),
#' generalized boosted models (`gbm`), and generalized additive models (`gam`).
#' Optionally, the function can perform cross-validated hyperparameter tuning via
#' the \pkg{caret} package when \code{tune = TRUE}.
#'
#' @param data A \code{data.frame} containing the treatment, covariate, and outcome variables.
#' @param treatment A single numeric column index of the treatment variable in \code{data}.
#' @param treatment_ref An optional reference level for the treatment factor.
#'   If not specified, the first level of the treatment factor is used as the reference.
#' @param covariate A numeric vector of column indices corresponding to the covariate variables.
#' @param gps_model A character string specifying the GPS estimation model to be used.
#'   Choices include:
#'   \itemize{
#'     \item \code{"logit"} – Multinomial logistic regression (via \pkg{nnet})
#'     \item \code{"gbm"} – Generalized boosted model (via \pkg{gbm})
#'     \item \code{"rf"} – Random forest (future extension)
#'     \item \code{"gam"} – Generalized additive model (via \pkg{VGAM})
#'   }
#' @param gps_params An optional named list of model-specific parameters.
#'   Each element corresponds to one model and is itself a list of hyperparameters.
#'   For example:
#'   \preformatted{
#'   list(
#'     logit = list(maxit = 200, decay = 1e-3, trace = FALSE),
#'     gam   = list(df_max = 6, maxit = 80, trace = FALSE),
#'     gbm   = list(n.trees = 4000, interaction.depth = 4,
#'                  shrinkage = 0.01, cv.folds = 5,
#'                  n.minobsinnode = 5, n.cores = 2)
#'   )
#'   }
#' @param tune Logical; if \code{TRUE}, performs hyperparameter tuning using
#'   cross-validation through \pkg{caret::train}.
#' @param tune_control Optional object created by \code{caret::trainControl()},
#'   specifying resampling and cross-validation options. If \code{NULL},
#'   a default 5-fold cross-validation setting is used.
#' @param tune_grids An optional named list of tuning grids, each defined as
#'   a \code{data.frame} whose column names correspond to valid tuning parameters
#'   for the selected model. Use \code{caret::getModelInfo("<method>")} to view
#'   the allowed parameter names.
#' @param seed Integer; random seed for reproducibility.
#'
#' @return A \code{data.frame} containing:
#'   \itemize{
#'     \item The original input data.
#'     \item Estimated GPS columns, named \code{gps_<treatment_level>}.
#'     \item Log-transformed GPS columns (log-odds ratios relative to the
#'           reference level), named \code{loggps_<treatment_level>}.
#'   }
#'
#' @details
#' The function automatically clips predicted probabilities to \eqn{[10^{-6}, 1 - 10^{-6}]}
#' to avoid numerical instability when computing log-ratios. When
#' \code{tune = TRUE}, the function validates that the provided
#' \code{tune_grids} contain only model-specific tuning parameters supported by
#' \pkg{caret}. If invalid parameters are found, an informative message is returned
#' indicating the correct parameter set.
#'
#' @importFrom splines bs
#' @importFrom nnet multinom
#' @importFrom gbm gbm gbm.perf
#' @importFrom VGAM vglm multinomial
#' @importFrom caret train trainControl getModelInfo
#' @examples
#' \dontrun{
#' gps_pre_process(
#'   data = df,
#'   treatment = 1,
#'   covariate = 2:6,
#'   gps_model = "gbm",
#'   tune = TRUE,
#'   tune_grids = list(
#'     gbm = expand.grid(
#'       n.trees = c(1000, 3000),
#'       interaction.depth = c(2, 3),
#'       shrinkage = c(0.01, 0.005),
#'       n.minobsinnode = c(5, 10)
#'     )
#'   )
#' )
#' }
gps_pre_process <- function(data,
                            treatment,
                            treatment_ref = NULL,
                            covariate,
                            gps_model = c("logit","gbm","gam"),
                            gps_params = NULL,
                            tune = FALSE,
                            tune_control = NULL,
                            tune_grids = NULL,
                            seed = 12345) {

  ## -- helpers --
  .merge_defaults <- function(defaults, user) {
    if (is.null(user)) return(defaults)
    stopifnot(is.list(user)); for (nm in names(user)) defaults[[nm]] <- user[[nm]]; defaults
  }
  .filter_to_formals <- function(fun, args) {
    keep <- intersect(names(args), names(formals(fun))); dropped <- setdiff(names(args), keep)
    list(keep = args[keep], dropped = dropped)
  }
  get_user <- function(name) { if (is.null(gps_params) || !is.list(gps_params)) return(NULL); gps_params[[name]] }
  .default_trainControl <- function() {
    caret::trainControl(method = "cv", number = 5,
                        classProbs = TRUE,
                        summaryFunction = caret::multiClassSummary,
                        savePredictions = "final",
                        allowParallel = TRUE)
  }
  .default_grids <- function(model_name) {
    switch(model_name,
           "logit" = expand.grid(decay = 10^seq(-4,-1,length.out=5)),
           "gbm"   = expand.grid(n.trees=c(1500,3000,4500),
                                 interaction.depth=c(2,3,4),
                                 shrinkage=c(0.01,0.005),
                                 n.minobsinnode=c(5,10)),
           NULL)
  }
  .validate_tunegrid <- function(method, grid) {
    mi <- caret::getModelInfo(method)[[method]]
    expected <- mi$parameters$parameter
    if (!is.data.frame(grid) || !setequal(colnames(grid), expected)) {
      stop(sprintf('Invalid tuneGrid for method "%s"; please use getModelInfo("%s")$%s$parameters to get proper parameters for tuning',
                   method, method, method))
    }
    invisible(TRUE)
  }

  ## -- input checks --
  if (!is.numeric(treatment) || length(treatment) != 1) stop("`treatment` option must be single numeric column number")
  if (!is.numeric(covariate) || length(covariate) < 1) stop("`covariate` option must be numeric column numbers for more than one column")
  p <- ncol(data)
  if (any(c(treatment, covariate) < 1 | c(treatment, covariate) > p))
    stop("treatment or covariate column number is less or greater than the total column number")

  trt_var  <- names(data)[treatment]
  cov_vars <- names(data)[covariate]
  trt_fac <- factor(data[[trt_var]])
  if (!is.null(treatment_ref)) {
    if (!treatment_ref %in% levels(trt_fac)) stop("`treatment_ref` must be a level in treatment column")
    trt_fac <- stats::relevel(trt_fac, ref = treatment_ref)
  }
  data[[trt_var]] <- trt_fac
  K <- nlevels(data[[trt_var]])
  if (K < 2) stop("treatment column must have at least 2 levels")

  model <- match.arg(gps_model)
  set.seed(seed)

  ## -- logit (multinomial) --
  if (model == "logit") {
    form <-  stats::reformulate(cov_vars, response = trt_var)
    if (isTRUE(tune)) {
      tc   <- if (is.null(tune_control)) .default_trainControl() else tune_control
      grid <- if (!is.null(tune_grids) && !is.null(tune_grids$logit)) tune_grids$logit else .default_grids("logit")
      .validate_tunegrid("multinom", grid)
      user_par <- get_user("logit"); user_par <- if (is.null(user_par)) list() else user_par
      fit <- caret::train(
        form, data = data, method = "multinom",
        metric = "logLoss", trControl = tc, tuneGrid = grid,
        trace = FALSE,
        MaxNWts = user_par$MaxNWts %||% 10000,
        maxit   = user_par$maxit %||% 200
      )
      gps <- stats::predict(fit, newdata = data, type = "prob")
      gps <- as.matrix(gps)[, levels(trt_fac), drop = FALSE]
    } else {
      logit_def <- list(trace = FALSE)
      logit_par <- .merge_defaults(logit_def, get_user("logit"))
      filt <- .filter_to_formals(nnet::multinom, logit_par)
      fit  <- do.call(nnet::multinom, c(list(formula = form, data = data), filt$keep))
      gps  <- stats::fitted(fit)
      if (is.null(dim(gps)) || ncol(gps) == 1) gps <- cbind(1 - gps, gps)
      colnames(gps) <- levels(trt_fac)
    }
  }

  ## -- gbm (multinomial) --
  else if (model == "gbm") {
    form <-  stats::reformulate(cov_vars, response = trt_var)
    if (isTRUE(tune)) {
      tc   <- if (is.null(tune_control)) .default_trainControl() else tune_control
      grid <- if (!is.null(tune_grids) && !is.null(tune_grids$gbm)) tune_grids$gbm else .default_grids("gbm")
      .validate_tunegrid("gbm", grid)
      user_par <- get_user("gbm"); user_par <- if (is.null(user_par)) list() else user_par
      fit <- caret::train(
        form, data = data, method = "gbm",
        distribution = "multinomial",
        metric = "logLoss", trControl = tc, tuneGrid = grid, verbose = FALSE
      )
      gps <- stats::predict(fit, newdata = data, type = "prob")
      gps <- as.matrix(gps)[, levels(trt_fac), drop = FALSE]
    } else {
      gbm_def <- list(n.trees=3000L, interaction.depth=3L, shrinkage=0.01,
                      cv.folds=5L, n.minobsinnode=10L, n.cores=1L, verbose=FALSE)
      gbm_par <- .merge_defaults(gbm_def, get_user("gbm"))
      filt <- .filter_to_formals(gbm::gbm, gbm_par)
      fit <- suppressWarnings(do.call(gbm::gbm,
                                      c(list(formula = form, data = data, distribution = "multinomial"),
                                        filt$keep)))
      best_iter <- if (isTRUE(gbm_par$cv.folds > 1L)) gbm::gbm.perf(fit, method = "cv", plot.it = FALSE) else gbm_par$n.trees
      pred_arr <- gbm::predict.gbm(object = fit, newdata = data, n.trees = best_iter, type = "response")
      gps <- if (is.array(pred_arr) && length(dim(pred_arr)) == 3L) pred_arr[, , 1, drop = FALSE][, , 1] else as.matrix(pred_arr)
      colnames(gps) <- levels(trt_fac)
    }
  }

  ## -- gam (VGAM multinomial) --
  else if (model == "gam") {
    make_term <- function(v) { x <- data[[v]]; u <- length(unique(x)); if (u <= 3L) v else sprintf("splines::bs(%s, df=%d)", v, min(5L, u - 1L)) }
    rhs  <- vapply(cov_vars, make_term, character(1))
    form <- stats::as.formula(paste(trt_var, "~", paste(rhs, collapse = " + ")))
    fit  <- suppressWarnings(VGAM::vglm(formula = form, family = VGAM::multinomial(), data = data))
    gps  <- stats::predict(fit, type = "response")
    colnames(gps) <- levels(trt_fac)
  }

  ## -- finalize outputs --
  clip_prob <- function(p, eps = 1e-6) pmax(pmin(p, 1 - eps), eps)
  gps <- clip_prob(as.matrix(gps), eps = 1e-6)
  colnames(gps) <- paste0("gps_", levels(data[[trt_var]]))
  loggps <- log(gps[, -1, drop = FALSE] / gps[, 1])
  colnames(loggps) <- paste0("loggps_", levels(data[[trt_var]])[-1])
  cbind(data, gps, loggps)
}
