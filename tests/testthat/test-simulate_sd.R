test_that("simulate_function_sd basic functionality works", {
  df <- data.frame(
    Date = c(1, 2, 3, 4),
    Diam = c(30, 33, 35, 38),
    CT   = c(0, 1, 0, 1)
  )

  res <- simulate_function_sd(df, n_sim = 10)
  expect_type(res, "list")
  expect_named(res, c("mean", "sd", "n_valid"))
  expect_true(is.numeric(res$mean))
  expect_true(is.numeric(res$sd))
})

test_that("simulate_function_sd errors for invalid input types", {
  df_bad <- data.frame(Date = c(1, 2), Diam = c(10, 12))
  expect_error(simulate_function_sd(df_bad), "must contain columns")

  df <- data.frame(Date = c("a", "b"), Diam = c(10, 12), CT = c(0, 1))
  expect_error(simulate_function_sd(df), "must be numeric")

  df2 <- data.frame(Date = c(1, 2), Diam = c(10, 12), CT = c(0, 2))
  expect_error(simulate_function_sd(df2), "must contain only 0 and 1")
})

test_that("simulate_function_sd reproducibility with seed", {
  df <- data.frame(
    Date = 1:4,
    Diam = c(30, 31, 33, 36),
    CT   = c(0, 1, 0, 1)
  )
  res1 <- simulate_function_sd(df, seed = 123, n_sim = 10)
  res2 <- simulate_function_sd(df, seed = 123, n_sim = 10)
  expect_equal(res1, res2)
})

test_that("simulate_function_sd warns when NAs are produced", {
  df <- data.frame(Date = 1:3, Diam = c(10, 10, 10), CT = c(0, 0, 0))
  expect_warning(
    res <- simulate_function_sd(df, n_sim = 10, seed = 123),
    "produced NA"
  )
  expect_true(is.list(res))
})

test_that("simulate_function_sd ignores infinite or NA FUN results", {
  df <- data.frame(Date = 1:4, Diam = c(30, 32, 34, 35), CT = c(0, 0, 1, 1))
  bad_fun <- function(d, x) ifelse(mean(x) > 0, sample(c(Inf,2), size =1), 1)
  expect_warning(res <- simulate_function_sd(df, FUN = bad_fun, n_sim = 10, seed = 123))
  expect_true(is.numeric(res$mean))
})

test_that("simulate_function_sd checks FUN structure and output", {
  df <- data.frame(Date = 1:4, Diam = c(10, 12, 14, 16), CT = c(0, 1, 0, 1))

  # Wrong number of arguments
  f_bad <- function(x) mean(x)
  expect_error(simulate_function_sd(df, FUN = f_bad), "two arguments")

  # Correct: single numeric output
  f_good <- function(d, x) mean(x)
  res <- simulate_function_sd(df, FUN = f_good, n_sim = 10)
  expect_true(is.numeric(res$mean))

  # Wrong: returns vector
  f_vec <- function(d, x) c(1, 2)
  expect_error(simulate_function_sd(df, FUN = f_vec, n_sim = 10))
})

test_that("simulate_function_sd stops when all results are invalid", {
  df <- data.frame(Date = 1:4, Diam = c(10, 12, 14, 16), CT = c(0, 1, 0, 1))
  f_bad <- function(d, x) NA_real_
  expect_error(simulate_function_sd(df, FUN = f_bad, n_sim = 10),
               "All simulations failed")
})
