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
                       hist_path = NULL)
{
  gps_model     <- match.arg(gps_model)
  outcome_model <- match.arg(outcome_model)

  # gps_pre_process
  dat_processed <- gps_pre_process(data, treatment, treatment_ref, covariate, gps_model)

  # overlap histogram
  if (hist == TRUE){
    hist = gps_histogram(dat_processed)
    ggplot2::ggsave(hist,file=hist_path,width=8,height=8)
  }

  ## making contrast matrix
  trt_lev  <- levels(factor(data[[names(data)[treatment]]]))
  contrast <- build_contrast(trt_lev, ref = treatment_ref)

  ## K-fold split
  n <- nrow(dat_processed)
  id <- sample.int(n)
  folds_id <- split(id, rep(1:folds, length.out = n))

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
                   contrast = contrast,nboot)
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
