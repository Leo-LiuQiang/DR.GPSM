test_that("gps_matching returns correct structure with simple input", {
  set.seed(123)

  df <- data.frame(
    trt = factor(rep(c("A", "B"), each = 3)),
    y   = c(1, 2, 3, 4, 5, 6),
    loggps_B = runif(6)
  )

  contrast <- matrix(c(-1, 1), nrow = 1, byrow = TRUE)
  colnames(contrast) <- c("A", "B")
  rownames(contrast) <- "BvA"

  pred <- matrix(0, nrow = 2, ncol = 6)

  out <- gps_matching(df, treatment = 1, outcome = 2,
                      pred = pred,
                      contrast = contrast,
                      nboot = 5)

  expect_type(out, "list")
  expect_named(out, c("estimate", "ci_lower", "ci_upper"))

  expect_equal(length(out$estimate), nrow(contrast))
  expect_equal(names(out$estimate), rownames(contrast))
})

test_that("gps_matching errors if loggps column is missing", {
  df <- data.frame(
    trt = factor(rep(c("A", "B"), each = 3)),
    y   = c(1, 2, 3, 4, 5, 6)
  )

  contrast <- matrix(c(-1, 1), nrow = 1)
  colnames(contrast) <- c("A", "B")
  rownames(contrast) <- "BvA"

  pred <- matrix(0, nrow = 2, ncol = 6)

  expect_error(
    gps_matching(df, treatment = 1, outcome = 2,
                 pred = pred, contrast = contrast, nboot = 5),
    regexp = "Can't find loggps"
  )
})

test_that("gps_matching works with 3 treatment levels", {
  set.seed(456)

  df <- data.frame(
    trt = factor(rep(c("A", "B", "C"), each = 2)),
    y   = c(1,2,3,4,5,6),
    loggps_B = runif(6),
    loggps_C = runif(6)
  )

  contrast <- matrix(c(
    -1,  1,  0,
    -1,  0,  1,
    0, -1,  1
  ), nrow = 3, byrow = TRUE)

  colnames(contrast) <- c("A", "B", "C")
  rownames(contrast) <- c("BvA", "CvA", "CvB")

  pred <- matrix(0, nrow = 3, ncol = 6)

  out <- gps_matching(df, treatment = 1, outcome = 2,
                      pred = pred, contrast = contrast, nboot = 5)

  expect_named(out, c("estimate", "ci_lower", "ci_upper"))
  expect_equal(length(out$estimate), 3)
})
