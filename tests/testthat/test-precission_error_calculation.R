
library(testthat)

# Testing measure_precision function ------------------------------------------

tb <- data.frame(
  ID = c("id1", "id1", "id1", "id1", "id2", "id2", "id2", "id2"),
  sensor = c("s1", "s2", "s1", "s2", "s1", "s2", "s1", "s2"),
  measurement = c(5, 6, 25, 30, 9, 8, 89, 86))


# Testing ERROR and MESSAGES

test_that("Error when missing columns", {
  tb1 <- tb
  colnames(tb1) <- c("id", "sens", "measurement")

  expect_error(
    ImgOMIC::measure_precision(tb1),

    "The following required columns are missing: ID, sensor"

    )
})


test_that("Error when NA values in ID", {
  tb1 <- rbind(tb, c(NA, "s1", 8))

  expect_error(
    ImgOMIC::measure_precision(tb1),

    "Your dataframe have 1 NAs in ID column"

  )
})

test_that("Error when no duplicates", {
  tb1 <- tb
  tb1$sensor <- c("s1", "s2", "s3", "s4", "s1", "s2" ,"s3", "s4")

  expect_error(
    ImgOMIC::measure_precision(tb1),

    "There are no duplicates in your dataframe"

  )
})

test_that("Message when there are values with NO-duplicate", {

  tb1 <- data.frame(
    ID = c("id1", "id1", "id1", "id1", "id2", "id2", "id2", "id2", "id2"),
    sensor = c("s1", "s2", "s3", "s4", "s1", "s2" ,"s3", "s4", "s4"),
    measurement = c(5, 6, 25, 30, 9, 8, 89, 86, 90))

  expect_message(
    ImgOMIC::measure_precision(tb1),

    "7 measurements have been removed since they do not have duplicate values"

  )
})

test_that("Error when ALL samples have NA in sensor or measurement columns", {
  tb1 <- tb
  tb1$sensor <- c(NA, NA, NA, NA, NA, NA, NA, NA)

  expect_error(
    ImgOMIC::measure_precision(tb1),

    "All IDs have NA values in measurement OR sensor columns"

  )
})


# Testing correct data corrections

test_that("Samples with NA values in sensor or measurement are removed", {

  tb1 <- data.frame(
    ID = c("id1", "id1", "id1", "id1", "id2", "id2", "id2", "id2"),
    sensor = c("s1", "s2", "s1", NA, "s1", "s2", "s1", "s2"),
    measurement = c(5, 6, 25, 30, 9, 8, 89, NA))

  res <- ImgOMIC::measure_precision(tb1)
  expect_equal(nrow(res), 1)
})

test_that("Measurement column is interpreted correctly even if it is character type", {

  tb1 <- data.frame(
    ID = c("id1", "id1", "id1", "id1", "id2", "id2", "id2", "id2"),
    sensor = c("s1", "s2", "s1", NA, "s1", "s2", "s1", "s2"),
    measurement = c(5, 6, 25, 30, 9, 8, 89, "kk"))

  res <- ImgOMIC::measure_precision(tb1)
  expect_equal(nrow(res), 1)
})


# Testing plot generation

test_that("Plots are generated when img = TRUE", {
  expect_silent(measure_precision(tb, img = TRUE))
})

test_that("No plots are generated when img = FALSE", {

  # Capture plot history before calling the function
  plot_history_before <- grDevices::recordPlot()

  # Run the function with img = FALSE
  measure_precision(tb, img = FALSE)

  # Capture plot history after calling the function
  plot_history_after <- grDevices::recordPlot()

  # Compare the two plot histories to ensure no new plots were added
  expect_identical(plot_history_before, plot_history_after)
})

# Testing correct MAD and RMS calculation

test_that("RMS is correctly calculated with just duplicates", {
  tb <- data.frame(
    ID = c("id1", "id1", "id1", "id1", "id2", "id2", "id2", "id2"),
    sensor = c("s1", "s2", "s1", "s2", "s1", "s2", "s1", "s2"),
    measurement = c(5, 6, 25, 30, 9, 8, 89, 86))

  res <- measure_precision(tb, img = FALSE)

  val1 <- as.double(res$RMS[res$sensor == "s1"])
  val2 <- sqrt(1/2*(((5-25)/sqrt(2)/mean(c(5, 25)))^2+((9-89)/sqrt(2)/mean(c(9, 89)))^2))*100
  expect_equal(val1, val2)

})

test_that("MAD is correctly calculated with just duplicates", {
  tb <- data.frame(
    ID = c("id1", "id1", "id1", "id1", "id2", "id2", "id2", "id2"),
    sensor = c("s1", "s2", "s1", "s2", "s1", "s2", "s1", "s2"),
    measurement = c(5, 6, 25, 30, 9, 8, 89, 86))

  res <- measure_precision(tb, img = FALSE)

  val1 <- as.double(res$MAD[res$sensor == "s1"])
  val2 <- sqrt(pi/2)*1/2*(abs((5-25)/sqrt(2)/mean(c(5, 25))) + abs((9-89)/sqrt(2)/mean(c(9, 89))))*100
  expect_equal(val1, val2)


})
