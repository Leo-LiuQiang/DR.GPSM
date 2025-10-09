#' GPS- or Covariate-matching ATE estimator  (numeric-column interface)
#'
#' @param data        Data frame (within-fold subset)
#' @param treatment   Column index of the treatment variable (factor)
#' @param outcome     Column index of the outcome variable
#' @param pred        K × n cross-fitted residual matrix from the outcome model
#' @param contrast    C(K, 2) × K contrast matrix from build_contrast()
#' @param nboot       Number of bootstrap replications
#' @param match_on    'gps' (default) or 'covariates'
#' @param covariate   Numeric vector of covariate column indices (required if match_on='covariates')
#' @param cov_distance 'euclidean' (default) or 'mahalanobis' (only for match_on='covariates')
#' @param standardize Logical; standardize covariate features before matching (default TRUE)
#' @param ridge       Small ridge for covariance in Mahalanobis whitening (default 1e-8)
#' @param match_ratio Integer ≥ 1; number of matches per unit per target group (default 1)
#'
#' @return list(estimate, ci_lower, ci_upper)
#' @export
gps_matching <- function(data,
                         treatment,
                         outcome,
                         pred,
                         contrast,
                         nboot,
                         match_on = c("gps","covariates"),
                         covariate = NULL,
                         cov_distance = c("mahalanobis","euclidean"),
                         standardize = FALSE,
                         ridge = 1e-8,
                         match_ratio = 1L) {

  match_on     <- match.arg(match_on)
  cov_distance <- match.arg(cov_distance)
  stopifnot(is.numeric(match_ratio), length(match_ratio) == 1L, match_ratio >= 1)

  n <- nrow(data)
  stopifnot(is.factor(data[[treatment]]))
  tlevel <- levels(data[[treatment]])
  K <- length(tlevel)

  index <- vector("list", K)
  mat   <- vector("list", K)

  if (match_on == "gps") {
    loggps_cols <- integer(K - 1L)
    for (k in 2:K) {
      need <- paste0("loggps_", tlevel[k])
      idx  <- match(need, names(data))
      if (is.na(idx))
        stop(sprintf("Can't find loggps %s (ensure gps_pre_process ran).", need))
      loggps_cols[k - 1L] <- idx
    }
    for (i in seq_len(K)) {
      index[[i]] <- which(data[, treatment] == tlevel[i])
      mat[[i]]   <- as.matrix(data[index[[i]], loggps_cols, drop = FALSE])
    }
  } else {
    if (is.null(covariate) || length(covariate) < 1L)
      stop("When match_on='covariates', please provide `covariate` column indices.")

    X <- stats::model.matrix(~ . - 1, data = data[, covariate, drop = FALSE])
    keep <- which(apply(X, 2, function(v) stats::var(v) > 0))
    if (length(keep) == 0L) stop("All covariate columns have zero variance in this fold.")
    X <- X[, keep, drop = FALSE]

    if (isTRUE(standardize)) {
      X <- scale(X); X[is.na(X)] <- 0
    }
    if (cov_distance == "mahalanobis") {
      S <- stats::cov(X)
      diag(S) <- diag(S) + ridge * mean(diag(S))
      L <- tryCatch(chol(solve(S)), error = function(e) NULL)
      if (is.null(L)) {
        eig <- eigen(S, symmetric = TRUE)
        L <- eig$vectors %*% diag(1 / sqrt(pmax(eig$values, ridge))) %*% t(eig$vectors)
      }
      X <- X %*% L
    }
    for (i in seq_len(K)) {
      index[[i]] <- which(data[, treatment] == tlevel[i])
      if (length(index[[i]]) == 0L)
        stop(sprintf("No samples for treatment level '%s' in this fold.", tlevel[i]))
      mat[[i]] <- as.matrix(X[index[[i]], , drop = FALSE])
    }
  }

  M_idx <- vector("list", K)
  for (i in seq_len(K)) {
    M_idx[[i]] <- vector("list", K)
    for (j in seq_len(K)) {
      if (i == j) {
        M_idx[[i]][[j]] <- matrix(index[[i]], nrow = length(index[[i]]), ncol = 1)
      } else {
        k_eff <- min(match_ratio, nrow(mat[[j]]))
        nn <- RANN::nn2(mat[[j]], mat[[i]], k = k_eff,
                        treetype = "kd", searchtype = "priority")$nn.idx
        if (is.null(dim(nn))) nn <- matrix(nn, ncol = 1)
        M_idx[[i]][[j]] <- matrix(index[[j]][nn], nrow = length(index[[i]]), ncol = k_eff)
      }
    }
  }

  group_of_i   <- as.integer(data[[treatment]])
  pos_in_group <- integer(n)
  for (g in seq_len(K)) pos_in_group[index[[g]]] <- seq_along(index[[g]])

  tau <- matrix(0, nrow = nrow(contrast), ncol = n)

  for (ii in 1:nrow(contrast)) {
    neg <- which(contrast[ii, ] == -1)
    pos <- which(contrast[ii, ] ==  1)

    for (i in 1:n) {
      j <- group_of_i[i]
      r <- pos_in_group[i]

      if (contrast[ii, j] == -1) {
        k <- pos
        mk <- M_idx[[j]][[k]][r, ]
        wk <- 1 / length(mk)

        tau[ii, mk] <- tau[ii, mk] + wk * (data[mk, outcome] - pred[k, mk])
        tau[ii, i ] <- tau[ii, i ] -       (data[i , outcome] - pred[k, i ])

      } else if (contrast[ii, j] == 1) {
        k <- neg
        mk <- M_idx[[j]][[k]][r, ]
        wk <- 1 / length(mk)

        tau[ii, mk] <- tau[ii, mk] - wk * (data[mk, outcome] - pred[k, mk])
        tau[ii, i ] <- tau[ii, i ] +       (data[i , outcome] - pred[k, i ])

      } else { # contrast[ii, j] == 0
        k1 <- neg; k2 <- pos
        m1 <- M_idx[[j]][[k1]][r, ]; w1 <- 1 / length(m1)
        m2 <- M_idx[[j]][[k2]][r, ]; w2 <- 1 / length(m2)

        tau[ii, m2] <- tau[ii, m2] + w2 * (data[m2, outcome] - pred[k2, m2])
        tau[ii, m1] <- tau[ii, m1] - w1 * (data[m1, outcome] - pred[k1, m1])
        tau[ii, i ] <- tau[ii, i ] -        pred[k1, i]       +  pred[k2, i]
      }
    }
  }

  est <- rowMeans(tau)
  tau <- sweep(tau, 1, est, "-")

  btstrp <- matrix(nrow = nboot, ncol = nrow(contrast))
  for (b in 1:nboot) {
    w <- c(stats::rmultinom(1, 2*n, rep(1/n, n)))
    btstrp[b, ] <- rowSums(w * tau) / (2*n)
  }

  cilower <- -apply(btstrp, 2, function(y) stats::quantile(y, 0.975)) + est
  ciupper <- -apply(btstrp, 2, function(y) stats::quantile(y, 0.025)) + est

  list(
    estimate = stats::setNames(as.numeric(est),   rownames(contrast)),
    ci_lower = stats::setNames(as.numeric(cilower), rownames(contrast)),
    ci_upper = stats::setNames(as.numeric(ciupper), rownames(contrast))
  )
}
