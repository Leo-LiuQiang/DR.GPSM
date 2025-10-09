#' Generalized Propensity Score Model fitting
#'
#' @param data Data frame containing treatment, covariates, and outcome variables
#' @param treatment Column index of the treatment variable
#' @param treatment_ref Reference treatment level (optional; default is the first level)
#' @param covariate Numeric vector of covariate column indices
#' @param gps_model Choice of GPS model: 'logit', 'rf', 'gbm', or 'gam'
#' @param gps_params     Optional named list of hyperparameters per GPS model.
#'                       For example: list(
#'                          logit = list(maxit=200, decay=1e-3, trace=FALSE),
#'                          gam   = list(df_max=6, maxit=80, trace=FALSE),
#'                          gbm = list(n.trees=4000, interaction.depth=4, shrinkage=0.01,
#'                                cv.folds=5, n.minobsinnode=5, n.cores=2))
#'
#' @return A data frame containing the original data plus estimated GPS columns (gps_) and log-transformed GPS columns (loggps_)
#' @export
#' @importFrom splines bs
gps_pre_process <- function(data,
                            treatment,
                            treatment_ref = NULL,
                            covariate,
                            gps_model=c("logit","gbm","gam"),
                            gps_params = NULL
                            ){
  # Column verification
  if (!is.numeric(treatment) || length(treatment) != 1)
    stop("`treatment` option must be single numeric column number")
  if (!is.numeric(covariate) || length(covariate) < 1)
    stop("`covariate` option must be numeric column numbers for more than one column")

  # model parameter extraction
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

  get_user <- function(name) {
    if (is.null(gps_params) || !is.list(gps_params)) return(NULL)
    gps_params[[name]]
  }

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

    logit_def <- list(trace = FALSE)
    logit_par <- .merge_defaults(logit_def, get_user("logit"))

    filt <- .filter_to_formals(nnet::multinom, logit_par)

    fit <- do.call(nnet::multinom, c(list(formula = form, data = data), filt$keep))

    gps  <-  stats::fitted(fit)

    if (is.null(dim(gps)) || ncol(gps) == 1) gps <- cbind(1 - gps, gps)
  }

  # Generalized Boosted Models
  else if (gps_model == "gbm") {
    form <-  stats::reformulate(cov_vars, response = trt_var)

    gbm_def <- list(
      n.trees           = 3000L,
      interaction.depth = 3L,
      shrinkage         = 0.01,
      cv.folds          = 5L,
      n.minobsinnode    = 10L,
      n.cores           = 1L,
      verbose           = FALSE
    )
    gbm_par <- .merge_defaults(gbm_def, get_user("gbm"))

    filt <- .filter_to_formals(gbm::gbm, gbm_par)
    fit <- suppressWarnings(do.call(gbm::gbm,
                   c(list(formula = form,
                          data    = data,
                          distribution = "multinomial"),
                     filt$keep)))

    best_iter <- if (isTRUE(gbm_par$cv.folds > 1L)) {
      gbm::gbm.perf(fit, method = "cv", plot.it = FALSE)
    } else {
      gbm_par$n.trees
    }

    pred_arr <- gbm::predict.gbm(object  = fit,
                                 newdata = data,
                                 n.trees = best_iter,
                                 type = "response")

    gps <- if (is.array(pred_arr) && length(dim(pred_arr)) == 3L) {
      pred_arr[, , 1, drop = FALSE][, , 1]
    } else {
      as.matrix(pred_arr)
    }
    colnames(gps) <- levels(trt_fac)
  }

  # Generalized Additive Model
  else if (gps_model == "gam") {
    make_term <- function(v) {
      x <- data[[v]]
      u <- length(unique(x))
      if (u <= 3L) {
        return(v)
      } else {
        df_use <- min(5L, u - 1L)
        return(sprintf("splines::bs(%s, df=%d)", v, df_use))
      }
    }
    rhs <- vapply(cov_vars, make_term, character(1))
    form <- stats::as.formula(paste(trt_var, "~", paste(rhs, collapse = " + ")))

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
