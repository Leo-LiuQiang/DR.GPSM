#' Generalized Propensity Score Model fitting
#'
#' @param data Data frame containing treatment, covariates, and outcome variables
#' @param treatment Column index of the treatment variable
#' @param treatment_ref Reference treatment level (optional; default is the first level)
#' @param covariate Numeric vector of covariate column indices
#' @param gps_model Choice of GPS model: 'logit', 'rf', 'gbm', or 'gam'
#'
#' @return A data frame containing the original data plus estimated GPS columns (gps_) and log-transformed GPS columns (loggps_)
#' @export
gps_pre_process <- function(data,
                            treatment,
                            treatment_ref = NULL,
                            covariate,
                            gps_model=c("logit","rf","gbm","gam")){
  # Column verification
  if (!is.numeric(treatment) || length(treatment) != 1)
    stop("`treatment` option must be single numeric column number")
  if (!is.numeric(covariate) || length(covariate) < 1)
    stop("`covariate` option must be numeric column numbers for more than one column")

  p <- ncol(data)
  if (any(c(treatment, covariate) < 1 |
          c(treatment, covariate) > p))
    stop("treatment or covariate column number is less or greater than the total column number")

  # Column names extraction
  trt_var  <- names(data)[treatment]
  cov_vars <- names(data)[covariate]

  # Factorize and set reference treatement group
  trt_fac <- factor(data[[trt_var]])
  if (!is.null(treatment_ref)) {
    if (!treatment_ref %in% levels(trt_fac))
      stop("`treatment_ref` must be a level in treatment column")
    trt_fac <- stats::relevel(trt_fac, ref = treatment_ref)
  }
  data[[trt_var]] <- trt_fac
  K <- nlevels(data[[trt_var]])

  # Treatment verification
  if (K < 2) stop("treatment column must have at least 2 levels")

  # Generalized Propensity Score Model fitting
  model <- match.arg(gps_model)

  # Multinomial logistic model
  if (gps_model == "logit") {
    form <-  stats::reformulate(cov_vars, response = trt_var)

    fit  <- nnet::multinom(formula = form,
                           data = data,
                           trace = FALSE)

    gps  <-  stats::fitted(fit)
    if (is.null(dim(gps))) {
      gps <- cbind(1 - gps, gps)
    } else if (ncol(gps) == 1) {
      gps <- cbind(1 - gps, gps)
    }
  }

  # Random forest
  else if (gps_model == "rf") {
    form <-  stats::reformulate(cov_vars, response = trt_var)

    fit <- randomForest::randomForest(formula = form,
                                      data = data,
                                      ntree = 500)

    gps <- stats::predict(fit, type = "prob")
  }

  # Generalized Boosted Models
  else if (gps_model == "gbm") {
    form <-  stats::reformulate(cov_vars, response = trt_var)

    fit <- suppressWarnings(gbm::gbm(formula = form,
                        data = data,
                        distribution = "multinomial",
                        n.trees = 3000,
                        cv.folds = 5,
                        interaction.depth = 3,
                        shrinkage = 0.01,
                        n.minobsinnode = 10,
                        n.cores = 1,
                        verbose = FALSE))

    best_iter <- gbm::gbm.perf(fit, method = "cv", plot.it = FALSE)

    pred_arr <- gbm::predict.gbm(object  = fit,
                                 newdata = data,
                                 n.trees = best_iter,
                                 type = "response")

    gps <- matrix(pred_arr[,,1], ncol = K, dimnames = list(NULL, levels(trt_fac)))
  }

  # Generalized Additive Model
  else if (gps_model == "gam") {
    smooth_terms <- paste0("s(",cov_vars,")", collapse = "+")
    form <- stats::as.formula(
      paste(trt_var, "~", smooth_terms),
      env = asNamespace("VGAM"))

    fit <- suppressWarnings(VGAM::vglm(formula = form,
                      family  = VGAM::multinomial(),
                      data = data))

    gps <- stats::predict(fit, type = "response")
  }

  clip_prob <- function(p, eps = 1e-6) pmax(pmin(p, 1 - eps), eps)
  gps <- clip_prob(gps, eps = 1e-6)
  colnames(gps) <- paste0("gps_", levels(data[[trt_var]]))

  # loggps modification
  loggps <- log(gps[, -1, drop = FALSE] / gps[, 1])
  colnames(loggps) <- paste0("loggps_", levels(data[[trt_var]])[-1])

  # output new data
  out <- cbind(data, gps, loggps)
  return(out)
}
