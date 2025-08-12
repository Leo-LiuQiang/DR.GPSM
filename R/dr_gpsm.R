#' Doubly-robust GPS matching estimator
#'
#' @param data           Data frame with treatment, covariates, and outcome variables
#' @param treatment      Column index of the treatment variable
#' @param treatment_ref  Reference treatment level (default = first level)
#' @param covariate      Numeric vector of covariate column indices
#' @param outcome        Column index of the outcome variable
#' @param gps_model      Choice of GPS model: 'logit', 'rf', 'gbm', or 'gam'
#' @param outcome_model  Choice of outcome model for doubly-robust adjustment: 'none', 'lm', or 'rf'
#' @param folds          Number of folds for cross-fitting (default = 2)
#' @param nboot          Number of bootstrap replications (default = 500)
#' @param hist           Logical; TRUE to save GPS histogram
#' @param hist_path      Character string specifying file path (pdf, png, jpg) to save the histogram (used if hist = TRUE)
#' @param gps_params     Optional named list of hyperparameters per GPS model.
#'                       For example: list(
#'                          rf  = list(ntree=1000, mtry=3, sampsize=..., nodesize=2, maxnodes=..., nPerm=5),
#'                          gbm = list(n.trees=4000, interaction.depth=4, shrinkage=0.01,
#'                                cv.folds=5, n.minobsinnode=5, n.cores=2),
#'                          logit = list(maxit=200, decay=1e-3, trace=FALSE),
#'                          gam   = list(df_max=6, maxit=80, trace=FALSE))
#'
#' @return A list with components:
#' \describe{
#'   \item{estimate}{Named numeric vector of doubly-robust ATE estimates for each contrast}
#'   \item{ci_lower}{Named numeric vector of lower 95% confidence bounds}
#'   \item{ci_upper}{Named numeric vector of upper 95% confidence bounds}
#' }
#' @export
dr_gpsm <- function(data,
                       treatment,
                       treatment_ref = NULL,
                       covariate,
                       outcome,
                       gps_model = c("logit","rf","gbm","gam"),
                       outcome_model = c("none","lm","rf"),
                       folds = 2,
                       nboot = 500,
                       hist = FALSE,
                       hist_path = NULL,
                       gps_params = NULL
                    )
{
  gps_model     <- match.arg(gps_model)
  outcome_model <- match.arg(outcome_model)

  # gps_pre_process
  dat_processed <- gps_pre_process(data, treatment, treatment_ref, covariate, gps_model, gps_params = gps_params)

  # overlap histogram
  if (isTRUE(hist)) {
    if (is.null(hist_path)) stop("When hist=TRUE, please provide hist_path.")
    p <- gps_histogram(dat_processed)
    ggplot2::ggsave(filename = hist_path, plot = p, width = 8, height = 8)
  }

  ## making contrast matrix
  trt_var  <- names(data)[treatment]
  trt_lev  <- levels(dat_processed[[trt_var]])
  contrast <- build_contrast(trt_lev, ref = treatment_ref)

  ## stratified K-fold split by treatment levels
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

  ## Outcome matrixes
  est_mat <- cil_mat <- ciu_mat <- matrix(NA_real_,
                                          nrow = nrow(contrast),
                                          ncol = folds)
  cov_vars <- names(data)[covariate]
  outcome_var <- names(data)[outcome]

  ## Loop for folds
  for (f in seq_along(folds_id)) {
    test  <- folds_id[[f]]
    train <- setdiff(seq_len(n), test)
    ## fit outcome models
    pred_test <- matrix(0,
                        nrow = ncol(contrast),
                        ncol = n)
    if (outcome_model != "none") {
      fit_fun <- switch(outcome_model,
                        rf = function(df) randomForest::randomForest(
                          stats::reformulate(cov_vars,
                                             response = outcome_var),
                          data = df),
                        lm = function(df) stats::lm(
                          stats::reformulate(cov_vars,
                                             response = outcome_var),
                          data = df))
      for (k in seq_along(trt_lev)) {
        idx_tr <- train[dat_processed[train, treatment] == trt_lev[k]]
        if (length(idx_tr) == 0L) next
        fit_k  <- fit_fun(dat_processed[idx_tr, , drop = FALSE])
        pred_test[k, test] <-
          stats::predict(fit_k, newdata = dat_processed[test, covariate, drop = FALSE])
      }
    }
    ## dr_gpsm
    g <- gps_matching(dat_processed[test, , drop = FALSE],
                   treatment, outcome,
                   pred = pred_test[, test, drop = FALSE],
                   contrast = contrast,nboot = nboot)
    est_mat[, f] <- g$estimate
    cil_mat[, f] <- g$ci_lower
    ciu_mat[, f] <- g$ci_upper
  }

  ## Final output
  out <- list(
    estimate = stats::setNames(rowMeans(est_mat),rownames(contrast)),
    ci_lower = stats::setNames(rowMeans(cil_mat),rownames(contrast)),
    ci_upper = stats::setNames(rowMeans(ciu_mat),rownames(contrast))
  )
  return(out)
}
