test_that("build_contrast works with default reference (first level)", {
  levs <- c("A", "B", "C")
  mat <- build_contrast(levs)

  expect_true(is.matrix(mat))

  expect_equal(dim(mat), c(3, 3))

  expect_true(all(mat %in% c(-1, 0, 1)))

  expect_true(all(grepl("v", rownames(mat))))

  expect_equal(colnames(mat), levs)
})

test_that("build_contrast uses specified reference level correctly", {
  levs <- c("A", "B", "C")
  mat_default <- build_contrast(levs)
  mat_Bref    <- build_contrast(levs, ref="B")

  expect_false(all(mat_default == mat_Bref))

  expect_equal(dim(mat_Bref), c(3, 3))
})

test_that("build_contrast errors if ref not in levels", {
  levs <- c("A", "B", "C")

  expect_error(build_contrast(levs, ref="Z"))
})

test_that("build_contrast works with 2 levels", {
  levs <- c("T1", "T2")
  mat <- build_contrast(levs)

  expect_equal(dim(mat), c(1, 2))
  expect_true(all(mat %in% c(-1, 0, 1)))
})
