test_that("gps_pre_process returns gps and loggps columns for logit model", {
  set.seed(123)
  df <- data.frame(
    trt = factor(rep(c("A","B"), each=5)),
    x1  = rnorm(10),
    x2  = runif(10)
  )

  out <- gps_pre_process(df, treatment=1, covariate=2:3, gps_model="logit")

  expect_true(any(grepl("^gps_", names(out))))
  expect_true(any(grepl("^loggps_", names(out))))

  gps_cols <- grep("^gps_", names(out), value = TRUE)
  expect_true(all(abs(rowSums(out[gps_cols]) - 1) < 1e-8))
})

test_that("gps_pre_process errors on invalid treatment or covariate input", {
  df <- data.frame(
    trt = factor(rep(c("A","B"), each=5)),
    x1  = rnorm(10),
    x2  = runif(10)
  )

  expect_error(gps_pre_process(df, treatment="trt", covariate=2:3, gps_model="logit"))

  expect_error(gps_pre_process(df, treatment=1, covariate=numeric(0), gps_model="logit"))

  expect_error(gps_pre_process(df, treatment=1, covariate=2:3, treatment_ref="Z", gps_model="logit"))

  df$trt <- factor(rep("A", 10))
  expect_error(gps_pre_process(df, treatment=1, covariate=2:3, gps_model="logit"))
})

test_that("gps_pre_process works with rf model if randomForest installed", {
  skip_if_not_installed("randomForest")

  set.seed(123)
  df <- data.frame(
    trt = factor(rep(c("A","B"), each=30)),
    x1  = rnorm(60),
    x2  = runif(60)
  )

  out <- gps_pre_process(df, treatment=1, covariate=2:3, gps_model="rf")
  expect_true(any(grepl("^gps_", names(out))))
})

test_that("gps_pre_process works with gbm model if gbm installed", {
  skip_if_not_installed("gbm")

  set.seed(123)
  df <- data.frame(
    trt = factor(rep(c("A","B"), each=30)),
    x1  = rnorm(60),
    x2  = runif(60)
  )

  out <- gps_pre_process(df, treatment=1, covariate=2:3, gps_model="gbm")
  expect_true(any(grepl("^gps_", names(out))))
})

test_that("gps_pre_process works with gam model if VGAM installed", {
  skip_if_not_installed("VGAM")

  set.seed(123)
  df <- data.frame(
    trt = factor(rep(c("A","B"), each=30)),
    x1  = rnorm(60),
    x2  = runif(60)
  )

  out <- gps_pre_process(df, treatment=1, covariate=2:3, gps_model="gam")
  expect_true(any(grepl("^gps_", names(out))))
})
