
library(testthat)

# Test mean_pairwise_diff -----------------------------------------------------

test_that("That mean_pairwise_diff have 2 outputs", {
  v <- c(52, 20)

  expect_length(mean_pairwise_diff(v), 2)

})

test_that("mean_pairwise_diff handles NA values", {
  expect_error(mean_pairwise_diff(c(1, 2, NA)), "There are NA values in the mean_pairwise_diff calculous")
})

test_that("mean_pairwise_diff handles single element", {
  result <- mean_pairwise_diff(c(5))
  expect_equal(result[1], NA_real_)  # First value should be NA
  expect_equal(result[2], 1)        # Length should be 1
})

test_that("mean_pairwise_diff handles two elements", {
  result <- mean_pairwise_diff(c(3, 6))
  expect_equal(result[1], (3 - 6) / sqrt(2))  # Pairwise difference
  expect_equal(result[2], 2)                    # Length should be 2
})

test_that("mean_pairwise_diff calculates correctly for multiple elements", {
  values <- c(1, 3, 6)
  result <- mean_pairwise_diff(values)
  combn_diffs <- combn(values, 2, function(x) abs(x[1] - x[2]) / sqrt(2))
  expect_equal(result[1], mean(combn_diffs))  # Mean of pairwise differences
  expect_equal(result[2], length(values))    # Length of input vector
})

# Test RMS_precision -----------------------------------------------------------

test_that("RMS_precision checks for missing columns", {
  df <- data.frame(mn = c(2, 5, 6))
  expect_error(RMS_precision(df), "The following required columns are missing: diff")

  df <- data.frame(diff = c(0.1, 0.5, 0.2))
  expect_error(RMS_precision(df), "The following required columns are missing: mn")
})

test_that("RMS_precision calculates correctly", {
  df <- data.frame(mn = c(2, 5, 6), diff = c(0.1, 0.5, 0.2))
  expected_RMS <- sqrt(sum((df$diff / df$mn)^2) / nrow(df)) * 100
  expect_equal(RMS_precision(df), expected_RMS)
})

# Test MAD_precision -----------------------------------------------------------

test_that("MAD_precision checks for missing columns", {
  df <- data.frame(mn = c(2, 5, 6))
  expect_error(MAD_precision(df), "The following required columns are missing: diff")

  df <- data.frame(diff = c(0.1, 0.5, 0.2))
  expect_error(MAD_precision(df), "The following required columns are missing: mn")
})

test_that("MAD_precision calculates correctly", {
  df <- data.frame(mn = c(2, 5, 6), diff = c(0.1, 0.5, 0.2))
  expected_MAD <- sqrt(pi / 2) * sum(abs(df$diff / df$mn)) / nrow(df) * 100
  expect_equal(MAD_precision(df), expected_MAD)
})
