test_that("gps_histogram returns a ggplot object with valid input", {
  df <- data.frame(
    gps_A = runif(50, 0.1, 0.9),
    gps_B = runif(50, 0.1, 0.9)
  )

  p <- gps_histogram(df)

  expect_s3_class(p, "ggplot")
})

test_that("gps_histogram errors when less than 2 gps columns are present", {
  df <- data.frame(
    gps_A = runif(10)
  )
  expect_error(gps_histogram(df))
})

test_that("gps_histogram respects custom palette", {
  df <- data.frame(
    gps_A = runif(20, 0.2, 0.8),
    gps_B = runif(20, 0.2, 0.8),
    gps_C = runif(20, 0.2, 0.8)
  )

  custom_palette <- c("red", "blue", "green")

  p <- gps_histogram(df, palette = custom_palette)

  fill_scale <- p$scales$get_scales("fill")
  expect_equal(fill_scale$palette(3), custom_palette)
})
