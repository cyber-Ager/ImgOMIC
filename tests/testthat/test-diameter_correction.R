# tests/testthat/test-correct-diameters.R

library(testthat)
library(ImgOMIC)  # replace with your actual package name


# Sample data for single patient ----------------------------------------------
data_multi <- data.frame(
  ID = as.factor(c("1", "1", "1", "1", "1", "1",  "2", "2", "2", "2")),
  Date = c(2015, 2016, 2017, 2020, 2021, 2024, 2013, 2016, 2018, 2019),
  Diam = c(30, 35, 34, 53, 50, 52, 38, 42, 50, 53),
  CT = as.factor(c(1, 0, 0, 1, 0, 0, 0, 0, 0, 1))
)

df_patient <- subset(data_multi, ID == "1")

# Test that correct_diameter_single returns expected structure
test_that("correct_diameter_single returns list with correct elements", {
  res <- correct_diameter_single(df_patient, sdUS = 1, sdCT = 1, sp = 2,
                                 dlim_inf = 0, dlim_sup = 100)
  expect_type(res, "list")
  expect_named(res, c("curves", "max_prob", "num_curves"))
  expect_true(is.numeric(res$max_prob) && res$max_prob >= 0)
  expect_true(is.integer(res$num_curves) || is.numeric(res$num_curves))
  expect_true(is.list(res$curves))
})

# TEST VALUES
test_that("correct_diameter_single returns expected curves and probability", {

  expected_curves <- list(
    c(30, 32, 34, 50, 50, 52),
    c(30, 35, 37, 50, 50, 52)
  )

  expected_max_prob <- 0.0001081654
  expected_num_curves <- 2

  # Check number of curves
  expect_equal(length(res$curves), expected_num_curves)

  # Check individual curves (order-insensitive)
  for (ec in expected_curves) {
    expect_true(any(sapply(res$curves, function(rc) all(abs(rc - ec) < 1e-6))))
  }

  # Check max_prob and num_curves
  expect_equal(res$max_prob, expected_max_prob, tolerance = 1e-5)
  expect_equal(res$num_curves, expected_num_curves)
})

# TEST ERRORS
test_that("correct_diameter_single throws error for empty data", {
  expect_error(
    correct_diameter_single(data.frame()),
    "Input data frame is empty"
  )
})

test_that("correct_diameter_single throws error for missing required columns", {
  df_invalid <- data.frame(Date = c(2020, 2021), Diam = c(30, 32))
  expect_error(
    correct_diameter_single(df_invalid),
    "Input data frame is missing required columns: CT"
  )
})




# Sample data for multiple patients---------------------------------------------


# Test that correct_diameters_all returns a data frame with a 'corrected' column

test_that("correct_diameters_all adjusts Diam and adds corrected flag", {
  res_all <- correct_diameters_all(data_multi, sdUS = 1, sdCT = 1, sp = 2,
                                   dlim_inf = 0, dlim_sup = 100)
  expect_s3_class(res_all, "data.frame")
  expect_true("corrected" %in% names(res_all))
  expect_equal(nrow(res_all), nrow(data_multi))
  # Rows with NA Diam should remain NA and corrected==FALSE
  na_rows <- is.na(data_multi$Diam)
  expect_true(all(is.na(res_all$Diam[na_rows]) & !res_all$corrected[na_rows]))
  # For patient A with two valid measures, corrected should be TRUE
  a_rows <- which(data_multi$ID == "A")
  expect_true(all(res_all$corrected[a_rows]))
})

# TEST PROCESSES or NOT
test_that("corrected flag is FALSE when <2 valid diameters or incomplete data", {
  # Case 1: Only one valid diameter
  df1 <- data.frame(
    ID = factor("A"),
    Date = 2015,
    Diam = 30,
    CT   = factor(1)
  )
  res1 <- correct_diameters_all(df1)
  expect_false(res1$corrected)

  # Case 2: Three entries but CT and Date have missing values
  df2 <- data.frame(
    ID = factor(rep("B", 3)),
    Date = c(2015, NA, 2017),
    Diam = c(30, 34, 36),
    CT   = factor(c(1, 0, NA))
  )
  res2 <- correct_diameters_all(df2)
  expect_true(all(!res2$corrected))

  # Case 3: Three complete valid diameters
  df3 <- data.frame(
    ID = factor(rep("C", 3)),
    Date = c(2015, 2016, 2017),
    Diam = c(30, 34, 36),
    CT   = factor(c(1, 0, 0))
  )
  res3 <- correct_diameters_all(df3)
  expect_true(all(res3$corrected))
})

# TEST average Curve
test_that("correct_diameters_all returns average of curves from correct_diameter_single", {

  data_multi <- data.frame(
    ID = as.factor(c("1", "1", "1", "1", "1", "1", "2", "2", "2", "2")),
    Date = c(2015, 2016, 2017, 2020, 2021, 2024, 2013, 2016, 2018, 2019),
    Diam = c(30, 35, 34, 53, 50, 52, 38, 42, 50, 53),
    CT = as.factor(c(1, 0, 0, 1, 0, 0, 0, 0, 0, 1))
  )

  df_patient <- subset(data_multi, ID == "1")

  # Get single-patient result
  res_single <- correct_diameter_single(
    df_patient,
    sdUS = 1, sdCT = 1, sp = 2,
    dlim_inf = 0, dlim_sup = 100
  )

  # Compute average curve
  expected_avg <- colMeans(do.call(rbind, res_single$curves))

  # Get result from correct_diameters_all
  res_all <- correct_diameters_all(
    df_patient,
    sdUS = 1, sdCT = 1, sp = 2,
    dlim_inf = 0, dlim_sup = 100
  )

  # Extract only the corrected diameters
  corrected_diam <- res_all$Diam[res_all$corrected]

  # Expect equal with tolerance due to floating point arithmetic
  expect_equal(corrected_diam, expected_avg, tolerance = 1e-6)
})


# TEST ERRORS
test_that("correct_diameters_all throws error for empty data", {
  expect_error(
    correct_diameters_all(data.frame()),
    "`data` must be a non-empty data frame."
  )
})

test_that("correct_diameters_all throws error for missing columns", {
  df_missing <- data.frame(ID = 1:3, Date = c(2020, 2021, 2022))
  expect_error(
    correct_diameters_all(df_missing),
    "Input data frame is missing required columns: Diam, CT"
  )
})

