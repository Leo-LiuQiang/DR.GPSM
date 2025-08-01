test_that("dr_gpsm runs correctly with minimal settings", {
  set.seed(123)
  df <- data.frame(
    trt = factor(rep(c("A","B"), each = 3)),
    x1  = rnorm(6),
    x2  = runif(6),
    y   = rnorm(6)
  )

  out <- dr_gpsm(
    data = df,
    treatment = 1,
    treatment_ref = "A",
    covariate = 2:3,
    outcome = 4,
    gps_model = "logit",
    outcome_model = "none",
    folds = 2,
    nboot = 5
  )

  expect_type(out, "list")
  expect_named(out, c("estimate", "ci_lower", "ci_upper"))

  expect_true(length(out$estimate) == 1)
})

test_that("dr_gpsm can generate histogram when hist=TRUE", {
  set.seed(123)
  df <- data.frame(
    trt = factor(rep(c("A","B"), each = 3)),
    x1  = rnorm(6),
    x2  = runif(6),
    y   = rnorm(6)
  )

  tmpfile <- tempfile(fileext = ".png")
  dr_gpsm(
    data = df,
    treatment = 1,
    treatment_ref = "A",
    covariate = 2:3,
    outcome = 4,
    gps_model = "logit",
    outcome_model = "none",
    folds = 2,
    nboot = 5,
    hist = TRUE,
    hist_path = tmpfile
  )

  expect_true(file.exists(tmpfile))
})

test_that("dr_gpsm errors on invalid gps_model or outcome_model", {
  set.seed(123)
  df <- data.frame(
    trt = factor(rep(c("A","B"), each = 3)),
    x1  = rnorm(6),
    x2  = runif(6),
    y   = rnorm(6)
  )

  expect_error(
    dr_gpsm(df, treatment=1, covariate=2:3, outcome=4,
            gps_model="invalid", outcome_model="none")
  )

  expect_error(
    dr_gpsm(df, treatment=1, covariate=2:3, outcome=4,
            gps_model="logit", outcome_model="invalid")
  )
})
