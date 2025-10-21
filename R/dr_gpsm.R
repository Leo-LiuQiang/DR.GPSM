#' Doubly-robust GPS matching estimator
#'
#' Implements a doubly-robust generalized propensity score (GPS) matching estimator
#' with optional cross-fitting and bootstrap inference. GPS can be estimated via
#' multinomial logistic regression, GBM, or VGAM-based GAM; outcome models can be
#' fitted via LM, RF, GBM, or GAM. Both GPS and outcome models support optional
#' caret-based hyperparameter tuning with parameter-grid validation.
#'
#' @param data Data frame with treatment, covariates, and outcome variables
#' @param treatment Column index of the treatment variable
#' @param treatment_ref Reference treatment level (default = first level)
#' @param covariate Numeric vector of covariate column indices
#' @param outcome Column index of the outcome variable
#' @param gps_model Choice of GPS model: 'logit', 'gbm', or 'gam'
#' @param outcome_model Choice of outcome model: 'none','lm','rf','gbm','gam'
#' @param folds Number of folds for cross-fitting (default = 2)
#' @param nboot Number of bootstrap replications (default = 500)
#' @param hist Logical; TRUE to save GPS histogram
#' @param hist_path File path (pdf/png/jpg) to save histogram if hist = TRUE
#' @param gps_params Optional named list of hyperparameters per GPS model
#' @param outcome_params Optional named list of hyperparameters per outcome model
#' @param match_on 'gps' (default) or 'covariates'
#' @param cov_distance For covariate matching: 'euclidean' or 'mahalanobis' (default 'mahalanobis')
#' @param standardize If TRUE, center/scale covariate features before matching (default FALSE)
#' @param match_ratio Integer >= 1; number of matches per unit per target group (default 1)
#' @param gps_tune Logical; if TRUE, enable caret tuning for GPS ('logit','gbm')
#' @param gps_tune_control Optional \code{caret::trainControl()} for GPS tuning; default is 5-fold CV
#' @param gps_tune_grids Optional named list of GPS tuning grids (e.g., \code{list(gbm = ..., logit = ...)})
#' @param gps_seed Integer seed for GPS tuning and fitting (default 2025)
#' @param outcome_tune Logical; if TRUE, enable caret tuning for outcome models ('rf','gbm')
#' @param outcome_tune_control Optional \code{caret::trainControl()} for outcome tuning
#' @param outcome_tune_grids Optional named list of outcome tuning grids (e.g., \code{list(rf = ..., gbm = ...)})
#' @param outcome_seed Integer seed for outcome tuning and fitting (default 2025)
#'
#' @return A list with components:
#' \describe{
#'   \item{estimate}{Named numeric vector of doubly-robust ATE estimates for each contrast.}
#'   \item{ci_lower}{Named numeric vector of lower 95% confidence bounds.}
#'   \item{ci_upper}{Named numeric vector of upper 95% confidence bounds.}
#' }
#'
#' @export
#' @importFrom stats na.omit
#' @importFrom mgcv gam
#' @importFrom caret train trainControl getModelInfo
dr_gpsm <- function(data,
                    treatment,
                    treatment_ref = NULL,
                    covariate,
                    outcome,
                    gps_model = c("logit","gbm","gam"),
                    outcome_model = c("none","lm","rf","gbm","gam"),
                    folds = 2,
                    nboot = 500,
                    hist = FALSE,
                    hist_path = NULL,
                    gps_params = NULL,
                    outcome_params = NULL,
                    match_on = c("gps","covariates"),
                    cov_distance = c("mahalanobis","euclidean"),
                    standardize = FALSE,
                    match_ratio = 1L,
                    gps_tune = FALSE,
                    gps_tune_control = NULL,
                    gps_tune_grids = NULL,
                    gps_seed = 2025,
                    outcome_tune = FALSE,
                    outcome_tune_control = NULL,
                    outcome_tune_grids = NULL,
                    outcome_seed = 2025) {

  ## -- helpers --
  .merge_defaults <- function(defaults, user) { if (is.null(user)) return(defaults); stopifnot(is.list(user)); for (nm in names(user)) defaults[[nm]] <- user[[nm]]; defaults }
  .filter_to_formals <- function(fun, args) { keep <- intersect(names(args), names(formals(fun))); list(keep = args[keep], dropped = setdiff(names(args), keep)) }
  get_outcome_user <- function(name) { if (is.null(outcome_params) || !is.list(outcome_params)) return(NULL); outcome_params[[name]] }
  .default_tc <- function() caret::trainControl(method="cv", number=5, classProbs=TRUE, summaryFunction=caret::multiClassSummary, savePredictions="final", allowParallel=TRUE)
  .validate_tunegrid <- function(method, grid) {
    mi <- caret::getModelInfo(method)[[method]]
    expected <- mi$parameters$parameter
    if (!is.data.frame(grid) || !setequal(colnames(grid), expected)) {
      stop(sprintf('Invalid tuneGrid for method "%s"; please use getModelInfo("%s")$%s$parameters to get proper parameters for tuning', method, method, method))
    }
    invisible(TRUE)
  }
  `%||%` <- function(a,b) if (is.null(a)) b else a

  ## -- args and preprocessing --
  gps_model     <- match.arg(gps_model)
  outcome_model <- match.arg(outcome_model)
  match_on      <- match.arg(match_on)
  cov_distance  <- match.arg(cov_distance)

  dat_processed <- gps_pre_process(
    data = data,
    treatment = treatment,
    treatment_ref = treatment_ref,
    covariate = covariate,
    gps_model = gps_model,
    gps_params = gps_params,
    tune = gps_tune,
    tune_control = gps_tune_control %||% .default_tc(),
    tune_grids = gps_tune_grids,
    seed = gps_seed
  )

  if (isTRUE(hist)) {
    if (is.null(hist_path)) stop("When hist=TRUE, please provide hist_path.")
    p <- gps_histogram(dat_processed)
    ggplot2::ggsave(filename = hist_path, plot = p, width = 8, height = 8)
  }

  trt_var  <- names(data)[treatment]
  trt_lev  <- levels(dat_processed[[trt_var]])
  contrast <- build_contrast(trt_lev, ref = treatment_ref)

  n <- nrow(dat_processed)
  make_folds <- function(fac, K) {
    stopifnot(is.factor(fac))
    inds  <- split(seq_along(fac), fac)
    folds <- vector("list", K)
    for (ids in inds) {
      ids <- sample(ids)
      sp  <- split(ids, rep(1:K, length.out = length(ids)))
      for (f in seq_len(K)) folds[[f]] <- c(folds[[f]], sp[[f]])
    }
    lapply(folds, sort)
  }
  trt_fac  <- dat_processed[[trt_var]]
  folds_id <- make_folds(trt_fac, folds)

  est_mat <- cil_mat <- ciu_mat <- matrix(NA_real_, nrow = nrow(contrast), ncol = folds)
  cov_vars    <- names(data)[covariate]
  outcome_var <- names(data)[outcome]

  y <- dat_processed[[outcome_var]]
  is_binary <- { if (is.factor(y)) nlevels(y) == 2L else { u <- unique(stats::na.omit(y)); length(u) == 2L && all(sort(u) %in% c(0,1)) } }

  ## -- cross-fitting loop --
  for (f in seq_along(folds_id)) {
    test  <- folds_id[[f]]
    train <- setdiff(seq_len(n), test)

    pred_test <- matrix(0, nrow = ncol(contrast), ncol = n)

    if (outcome_model != "none") {

      ## -- build outcome fitter (with optional caret tuning for rf/gbm) --
      set.seed(outcome_seed)

      fit_fun <-
        switch(outcome_model,

               lm = {
                 function(df) stats::lm(stats::reformulate(cov_vars, response = outcome_var), data = df)
               },

               rf = {
                 if (isTRUE(outcome_tune)) {
                   tc   <- outcome_tune_control %||% .default_tc()
                   grid <- outcome_tune_grids$rf
                   if (!is.null(grid)) .validate_tunegrid("rf", grid) else {
                     grid <- expand.grid(mtry = floor(sqrt(length(cov_vars))) + c(-1,0,1))
                     grid$mtry <- pmax(1, grid$mtry)
                     .validate_tunegrid("rf", grid)
                   }
                   function(df) {
                     yy <- df[[outcome_var]]
                     is_cls <- is_binary
                     if (is_cls) yy <- factor(yy, levels = c(0,1), labels = c("zero","one"))
                     form <- stats::reformulate(cov_vars, response = outcome_var)
                     df2 <- df
                     df2[[outcome_var]] <- yy
                     caret::train(
                       form, data = df2,
                       method = "rf",
                       trControl = tc,
                       tuneGrid  = grid
                     )
                   }
                 } else {
                   rf_def <- list(ntree = 500L, mtry = NULL, nodesize = 1L)
                   rf_par <- .merge_defaults(rf_def, get_outcome_user("rf"))
                   function(df) {
                     form <- stats::reformulate(cov_vars, response = outcome_var)
                     do.call(randomForest::randomForest, c(list(formula = form, data = df), rf_par))
                   }
                 }
               },

               gbm = {
                 if (isTRUE(outcome_tune)) {
                   tc   <- outcome_tune_control %||% .default_tc()
                   grid <- outcome_tune_grids$gbm
                   if (!is.null(grid)) .validate_tunegrid("gbm", grid) else {
                     grid <- expand.grid(
                       n.trees = c(1500, 3000),
                       interaction.depth = c(2,3,4),
                       shrinkage = c(0.01, 0.005),
                       n.minobsinnode = c(5,10)
                     )
                     .validate_tunegrid("gbm", grid)
                   }
                   function(df) {
                     yy <- df[[outcome_var]]
                     is_cls <- is_binary
                     if (is_cls) yy <- factor(yy, levels = c(0,1), labels = c("zero","one"))
                     form <- stats::reformulate(cov_vars, response = outcome_var)
                     df2 <- df
                     df2[[outcome_var]] <- yy
                     caret::train(
                       form, data = df2,
                       method = "gbm",
                       distribution = if (is_binary) "bernoulli" else "gaussian",
                       metric = if (is_binary) "logLoss" else "RMSE",
                       trControl = tc,
                       tuneGrid  = grid,
                       verbose   = FALSE
                     )
                   }
                 } else {
                   gbm_def <- list(n.trees=3000L, interaction.depth=3L, shrinkage=0.01,
                                   cv.folds=0L, n.minobsinnode=10L, n.cores=1L, verbose=FALSE)
                   gbm_par <- .merge_defaults(gbm_def, get_outcome_user("gbm"))
                   function(df) {
                     form <- stats::reformulate(cov_vars, response = outcome_var)
                     distn <- if (is_binary) "bernoulli" else "gaussian"
                     do.call(gbm::gbm, c(list(formula = form, data = df, distribution = distn), gbm_par))
                   }
                 }
               },

               gam = {
                 gam_def <- list(df_max = 6L, family = NULL, method = "REML")
                 gam_par <- .merge_defaults(gam_def, get_outcome_user("gam"))
                 rhs_terms <- vapply(cov_vars, function(v) {
                   x <- dat_processed[[v]]
                   if (is.numeric(x)) sprintf("s(%s, k=%d)", v, gam_par$df_max) else v
                 }, character(1))
                 fam_use <- if (!is.null(gam_par$family)) gam_par$family else if (is_binary) stats::binomial() else stats::gaussian()
                 function(df) {
                   form <- stats::as.formula(paste(outcome_var, "~", paste(rhs_terms, collapse = " + ")))
                   mgcv::gam(formula = form, family = fam_use, data = df, method = gam_par$method)
                 }
               }
        )

      ## -- fit within-treatment models and predict on test fold --
      for (k in seq_along(trt_lev)) {
        idx_tr <- train[dat_processed[train, treatment] == trt_lev[k]]
        if (length(idx_tr) == 0L) next
        df_k <- dat_processed[idx_tr, , drop = FALSE]
        fit_k <- fit_fun(df_k)

        if (outcome_model %in% c("rf","gbm") && isTRUE(outcome_tune)) {
          if (is_binary) {
            prob_mat <- stats::predict(fit_k, newdata = dat_processed[test, , drop = FALSE], type = "prob")
            p <- as.numeric(prob_mat[, "one"])
          } else {
            p <- as.numeric(stats::predict(fit_k, newdata = dat_processed[test, , drop = FALSE], type = "raw"))
          }
          pred_test[k, test] <- p
        } else if (outcome_model == "gbm") {
          ntrees <- if (!is.null(fit_k$cv.folds) && fit_k$cv.folds > 0L) gbm::gbm.perf(fit_k, method = "cv", plot.it = FALSE) else fit_k$n.trees
          p <- gbm::predict.gbm(fit_k, newdata = dat_processed[test, , drop = FALSE], n.trees = ntrees, type = "response")
          pred_test[k, test] <- as.numeric(p)
        } else if (outcome_model == "lm") {
          p <- stats::predict(fit_k, newdata = dat_processed[test, covariate, drop = FALSE], type = "response")
          pred_test[k, test] <- as.numeric(p)
        } else if (outcome_model == "rf") {
          p <- stats::predict(fit_k, newdata = dat_processed[test, covariate, drop = FALSE], type = "response")
          pred_test[k, test] <- as.numeric(p)
        } else if (outcome_model == "gam") {
          p <- stats::predict(fit_k, newdata = dat_processed[test, , drop = FALSE], type = "response")
          pred_test[k, test] <- as.numeric(p)
        }
      }
    }

    g <- gps_matching(dat_processed[test, , drop = FALSE],
                      treatment, outcome,
                      pred = pred_test[, test, drop = FALSE],
                      contrast = contrast, nboot = nboot,
                      match_on = match_on,
                      covariate = covariate,
                      cov_distance = cov_distance,
                      standardize = standardize,
                      match_ratio = match_ratio)

    est_mat[, f] <- g$estimate
    cil_mat[, f] <- g$ci_lower
    ciu_mat[, f] <- g$ci_upper
  }

  list(
    estimate = stats::setNames(rowMeans(est_mat), rownames(contrast)),
    ci_lower = stats::setNames(rowMeans(cil_mat), rownames(contrast)),
    ci_upper = stats::setNames(rowMeans(ciu_mat), rownames(contrast))
  )
}
