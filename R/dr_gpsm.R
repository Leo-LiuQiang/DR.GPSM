#' Doubly-robust GPS matching estimator
#'
#' @param data           Data frame with treatment, covariates, and outcome variables
#' @param treatment      Column index of the treatment variable
#' @param treatment_ref  Reference treatment level (default = first level)
#' @param covariate      Numeric vector of covariate column indices
#' @param outcome        Column index of the outcome variable
#' @param gps_model      Choice of GPS model: 'logit', 'gbm', or 'gam'
#' @param outcome_model  Choice of outcome model: 'none','lm','rf','gbm','gam'
#' @param folds          Number of folds for cross-fitting (default = 2)
#' @param nboot          Number of bootstrap replications (default = 500)
#' @param hist           Logical; TRUE to save GPS histogram
#' @param hist_path      File path (pdf/png/jpg) to save histogram if hist=TRUE
#' @param gps_params     Optional named list of hyperparameters per GPS model
#' @param outcome_params Optional named list of hyperparameters per outcome model
#'   For example: list(
#'     lm    = list(),  # (unused; here for symmetry)
#'     rf    = list(ntree = 500, mtry = NULL, nodesize = 1),
#'     gbm   = list(n.trees=3000, interaction.depth=3, shrinkage=0.01,
#'                  cv.folds=0, n.minobsinnode=10, n.cores=1, verbose=FALSE),
#'     gam   = list(df_max = 6, family = NULL, method = "REML")
#'   )
#' @param match_on     'gps' (default) or 'covariates'
#' @param cov_distance For covariate matching: 'euclidean' or 'mahalanobis' (default 'mahalanobis')
#' @param standardize  If TRUE, center/scale covariate features before matching (default FALSE)
#' @param match_ratio  Integer ≥ 1; number of matches per unit per target group (default 1)
#'
#' @return A list with components:
#' \describe{
#'   \item{estimate}{Named numeric vector of doubly-robust ATE estimates for each contrast}
#'   \item{ci_lower}{Named numeric vector of lower 95% confidence bounds}
#'   \item{ci_upper}{Named numeric vector of upper 95% confidence bounds}
#' }
#' @export
#' @importFrom stats na.omit
#' @importFrom mgcv gam
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
                    match_ratio = 1L) {

  .merge_defaults <- function(defaults, user) {
    if (is.null(user)) return(defaults)
    stopifnot(is.list(user))
    for (nm in names(user)) defaults[[nm]] <- user[[nm]]
    defaults
  }
  .filter_to_formals <- function(fun, args) {
    keep <- intersect(names(args), names(formals(fun)))
    dropped <- setdiff(names(args), keep)
    list(keep = args[keep], dropped = dropped)
  }
  get_outcome_user <- function(name) {
    if (is.null(outcome_params) || !is.list(outcome_params)) return(NULL)
    outcome_params[[name]]
  }

  gps_model     <- match.arg(gps_model)
  outcome_model <- match.arg(outcome_model)
  match_on      <- match.arg(match_on)
  cov_distance  <- match.arg(cov_distance)

  dat_processed <- gps_pre_process(data, treatment, treatment_ref, covariate,
                                   gps_model, gps_params = gps_params)

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
  is_binary <- {
    if (is.factor(y)) nlevels(y) == 2L else {
      u <- unique(na.omit(y))
      length(u) == 2L && all(sort(u) %in% c(0,1))
    }
  }

  for (f in seq_along(folds_id)) {
    test  <- folds_id[[f]]
    train <- setdiff(seq_len(n), test)

    pred_test <- matrix(0, nrow = ncol(contrast), ncol = n)
    if (outcome_model != "none") {

      fit_fun <-
        switch(outcome_model,

               lm = {
                 function(df) {
                   stats::lm(stats::reformulate(cov_vars, response = outcome_var), data = df)
                 }
               },

               rf = {
                 rf_def <- list(ntree = 500L, mtry = NULL, nodesize = 1L)
                 rf_par <- .merge_defaults(rf_def, get_outcome_user("rf"))
                 function(df) {
                   form <- stats::reformulate(cov_vars, response = outcome_var)
                   do.call(randomForest::randomForest,
                           c(list(formula = form, data = df), rf_par))
                 }
               },

               gbm = {
                 gbm_def <- list(
                   n.trees           = 3000L,
                   interaction.depth = 3L,
                   shrinkage         = 0.01,
                   cv.folds          = 0L,
                   n.minobsinnode    = 10L,
                   n.cores           = 1L,
                   verbose           = FALSE
                 )
                 gbm_par <- .merge_defaults(gbm_def, get_outcome_user("gbm"))
                 function(df) {
                   form <- stats::reformulate(cov_vars, response = outcome_var)
                   distn <- if (is_binary) "bernoulli" else "gaussian"
                   do.call(gbm::gbm,
                           c(list(formula = form, data = df, distribution = distn),
                             gbm_par))
                 }
               },

               gam = {
                 gam_def <- list(df_max = 6L, family = NULL, method = "REML")
                 gam_par <- .merge_defaults(gam_def, get_outcome_user("gam"))
                 rhs_terms <- vapply(cov_vars, function(v) {
                   x <- dat_processed[[v]]
                   if (is.numeric(x)) sprintf("s(%s, k=%d)", v, gam_par$df_max) else v
                 }, character(1))
                 fam_use <- if (!is.null(gam_par$family)) {
                   gam_par$family
                 } else {
                   if (is_binary) stats::binomial() else stats::gaussian()
                 }
                 function(df) {
                   form <- stats::as.formula(paste(outcome_var, "~", paste(rhs_terms, collapse = " + ")))
                   mgcv::gam(formula = form, family = fam_use, data = df, method = gam_par$method)
                 }
               }
        )

      for (k in seq_along(trt_lev)) {
        idx_tr <- train[dat_processed[train, treatment] == trt_lev[k]]
        if (length(idx_tr) == 0L) next
        fit_k  <- fit_fun(dat_processed[idx_tr, , drop = FALSE])

        if (outcome_model == "gbm") {
          ntrees <- if (!is.null(fit_k$cv.folds) && fit_k$cv.folds > 0L) {
            gbm::gbm.perf(fit_k, method = "cv", plot.it = FALSE)
          } else fit_k$n.trees
          p <- gbm::predict.gbm(fit_k,
                                newdata = dat_processed[test, , drop = FALSE],
                                n.trees = ntrees, type = "response")
          pred_test[k, test] <- as.numeric(p)
        } else {
          p <- stats::predict(fit_k, newdata = dat_processed[test, covariate, drop = FALSE], type = "response")
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
