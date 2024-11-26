
library(testthat)
library(dplyr)

# Testing measure_precision function ------------------------------------------

tb <- data.frame(
  ID = c("id1", "id1", "id1", "id1", "id2", "id2", "id2", "id2"),
  sensor = c("s1", "s2", "s1", "s2", "s1", "s2", "s1", "s2"),
  measurement = c(5, 6, 25, 30, 9, 8, 89, 86))


# Testing ERROR and MESSAGES

test_that("Error when missing columns in measure_precision", {
  tb1 <- tb
  colnames(tb1) <- c("id", "sens", "measurement")

  expect_error(
    ImgOMIC::measure_precision(tb1),

    "The following required columns are missing: ID, sensor"

    )
})


test_that("Error when NA values in ID column in measure_precision", {
  tb1 <- rbind(tb, c(NA, "s1", 8))

  expect_error(
    ImgOMIC::measure_precision(tb1),

    "Your dataframe have 1 NAs in ID column"

  )
})

test_that("Error when no duplicates values in measure_precision", {
  tb1 <- tb
  tb1$sensor <- c("s1", "s2", "s3", "s4", "s1", "s2" ,"s3", "s4")

  expect_error(
    ImgOMIC::measure_precision(tb1),

    "There are no duplicates in your dataframe"

  )
})

test_that("Message when in measure_precision there are values without duplicated value", {

  tb1 <- data.frame(
    ID = c("id1", "id1", "id1", "id1", "id2", "id2", "id2", "id2", "id2"),
    sensor = c("s1", "s2", "s3", "s4", "s1", "s2" ,"s3", "s4", "s4"),
    measurement = c(5, 6, 25, 30, 9, 8, 89, 86, 90))

  expect_message(
    ImgOMIC::measure_precision(tb1),

    "7 measurements have been removed since they do not have duplicate values"

  )
})

test_that("measure_precision | Error when ALL samples have NA in sensor or measurement columns", {
  tb1 <- tb
  tb1$sensor <- c(NA, NA, NA, NA, NA, NA, NA, NA)

  expect_error(
    ImgOMIC::measure_precision(tb1),

    "All IDs have NA values in measurement OR sensor columns"

  )
})


# Testing correct data corrections

test_that("measure_precision | Samples with NA values in sensor or measurement are removed", {

  tb1 <- data.frame(
    ID = c("id1", "id1", "id1", "id1", "id2", "id2", "id2", "id2"),
    sensor = c("s1", "s2", "s1", NA, "s1", "s2", "s1", "s2"),
    measurement = c(5, 6, 25, 30, 9, 8, 89, NA))

  res <- ImgOMIC::measure_precision(tb1)
  expect_equal(nrow(res), 1)
})

test_that("measure_precision | Measurement column is interpreted correctly even if it is character type", {

  tb1 <- data.frame(
    ID = c("id1", "id1", "id1", "id1", "id2", "id2", "id2", "id2"),
    sensor = c("s1", "s2", "s1", NA, "s1", "s2", "s1", "s2"),
    measurement = c(5, 6, 25, 30, 9, 8, 89, "kk"))

  res <- ImgOMIC::measure_precision(tb1)
  expect_equal(nrow(res), 1)
})


# Testing plot generation

test_that("measure_precision | Plots are generated when img = TRUE", {
  expect_silent(measure_precision(tb, img = TRUE))
})

test_that("measure_precision | No plots are generated when img = FALSE", {

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

test_that("measure_precision | RMS is correctly calculated with just duplicates", {
  tb <- data.frame(
    ID = c("id1", "id1", "id1", "id1", "id2", "id2", "id2", "id2"),
    sensor = c("s1", "s2", "s1", "s2", "s1", "s2", "s1", "s2"),
    measurement = c(5, 6, 25, 30, 9, 8, 89, 86))

  res <- measure_precision(tb, img = FALSE)

  val1 <- as.double(res$RMS[res$sensor == "s1"])
  val2 <- sqrt(1/2*(((5-25)/sqrt(2)/mean(c(5, 25)))^2+((9-89)/sqrt(2)/mean(c(9, 89)))^2))*100
  expect_equal(val1, val2)

})

test_that("measure_precision | MAD is correctly calculated with just duplicates", {
  tb <- data.frame(
    ID = c("id1", "id1", "id1", "id1", "id2", "id2", "id2", "id2"),
    sensor = c("s1", "s2", "s1", "s2", "s1", "s2", "s1", "s2"),
    measurement = c(5, 6, 25, 30, 9, 8, 89, 86))

  res <- measure_precision(tb, img = FALSE)

  val1 <- as.double(res$MAD[res$sensor == "s1"])
  val2 <- sqrt(pi/2)*1/2*(abs((5-25)/sqrt(2)/mean(c(5, 25))) + abs((9-89)/sqrt(2)/mean(c(9, 89))))*100
  expect_equal(val1, val2)


})

# Edge cases

test_that("measure_precision | Error when introducing empty dataframe in measuring precision", {
  tb1 <- data.frame()

  expect_error(
    measure_precision(tb1),

    "The following required columns are missing: ID, sensor, measurement"

  )

  tb1 <- data.frame(ID = numeric(0), sensor = numeric(0), measurement = numeric(0))

  expect_error(
    measure_precision(tb1),

    "All IDs have NA values in measurement OR sensor columns"

  )

})

test_that("measure_precision | Error introducing a data.frame with just one row in measuring precision", {
  tb1 <- data.frame(ID = "id1", sensor = "s1", measurement = 1)

  expect_error(
    measure_precision(tb1),

    "There are no duplicates in your dataframe"

  )

})

# Testing error_propagation function-------------------------------------------

df <- measure_precision(tb, img = FALSE)
colnames(df) <- c("sensor", "data", "RMS", "err", "n_precision")

fun <- expression(s1/s2)

# Testing ERROR or MESSAGES

test_that("Error when missing columns in error_propagation", {
  df1 <- df
  colnames(df1) <- c("X", "X", "RMS", "err", "n_precision")

  expect_error(
    error_propagation(df1, fun),

    "The following required columns are missing: sensor, data"

  )
})

test_that("error_propagation | Error when data column is not a df", {
  df1 <- df
  df1$data <- df1$n_precision

  expect_error(
    error_propagation(df1, fun),

    "The following required columns are missing: ID, mn, n"

  )
})

# For the next case I know that I am modifying just one of the data df-s, I mean I am getting just the one from the first row.
# In the future I can improve the code, but by the moment I dont think it would be necessary.
test_that("error_propagation | Error when data column is missing a required column", {
  df1 <- df
  colnames(df1$data[[1]]) <- c("X", "diff", "XX", "n")

  expect_error(
    error_propagation(df1, fun),

    "The following required columns are missing: ID, mn"

  )
})


test_that("error_propagation | Error when fun variable is missing", {

  expect_error(
    error_propagation(df),

    "fun is missing. It is neccessary to provide the fun expression."

  )
})

test_that("error_propagation | Error when fun variables are not shown in the input dataframe", {

  expect_error(
    error_propagation(df, expression(x + y)),

    "The following variables defined in the equation are not shown in the input tibble: x, y"

  )
})


test_that("error_propagation | Error when fun is not expresion type", {

  expect_error(
    error_propagation(df, c("hello", 6)),

    "fun need to be expression type"

  )
})

test_that("error_propagation | Tst that when NAs are introduced in the data df are being ignored", {

  df1 <- df
  df1$data[[1]] <- rbind(df1$data[[1]], c("id3", 0, NA, 0))

  expect_message(
    error_propagation(df1, fun),

    "1 samples have been excluded from the calulation due to missing values in sensor data."

  )

  res1 <- error_propagation(df, fun)
  res2 <- error_propagation(df1, fun)
  expect_equal(res1, res2)

})

test_that("error_propagation | Error when there is NO sample with ALL the reuired measurements", {
  df1 <- df
  df1$data[[1]]$mn <- c(NA, 49)
  df1$data[[2]]$mn <- c(18, NA)


  expect_error(
    error_propagation(df1, fun),

    "0 samples with ALL the needed measurements. Impossible to calculate the new measurement in any sample."

  )
})

test_that("error_propagation | Error when there is NA values in required sensor columns", {
  df1 <- rbind(df, df)
  df1$sensor <- c("s1", "s2", "s3", "s4")
  df1$err <- c(10, NA, NA, 5)

  fun1 <- expression(s1+s2+s3+s4)


  expect_error(
    error_propagation(df1, fun1),

    "The following sensors have NA values in err column: s2, s3"

  )
})

# Testing the correct CALCULATION of the parameters

test_that("error_propagation | Extra sensor rows are not affecting the final results", {
  df1 <- rbind(df, df)
  df1$sensor <- c("s1", "s2", "s3", "s4")

  res1 <- error_propagation(df1, fun)
  res2 <- error_propagation(df, fun)
  expect_equal(res1, res2)

})

test_that("error_propagation | Measure results are correct", {
  df1 <- rbind(df, df)
  df1$sensor <- c("s1", "s2", "s3", "s4")

  fun1 <- expression((s1*s2)/s3^2)

  res11 <- error_propagation(df1, fun1)$measure
  res12 <- (df1$data[[1]]$mn * df1$data[[2]]$mn)/df1$data[[3]]$mn^2

  expect_equal(res11, res12)

  fun2 <- expression(s1 + s4)

  res21 <- error_propagation(df1, fun2)$measure
  res22 <- (df1$data[[1]]$mn + df1$data[[4]]$mn)

  expect_equal(res21, res22)

})

test_that("error_propagation | Error results are correct", {

  fun1 <- expression((s1*s2))

  # dfun1/ds1 = s2 | dfun1/ds2 = s1

  res11 <- error_propagation(df, fun1)$error
  res12 <- sqrt((df$data[[2]]$mn * df$err[[1]]/100 * df$data[[1]]$mn / sqrt(2))^2 +
    (df$data[[1]]$mn * df$err[[2]]/100 * df$data[[2]]$mn / sqrt(2))^2)

  expect_equal(res11, res12)

  fun2 <- expression(s1 + 2*s2)

  # dfun2/ds1 = 1 | dfun2/ds2 = 2

  res21 <- error_propagation(df, fun2)$error
  res22 <- sqrt((1 * df$err[[1]]/100 * df$data[[1]]$mn / sqrt(2))^2 +
                  (2 * df$err[[2]]/100 * df$data[[2]]$mn / sqrt(2))^2)

  expect_equal(res21, res22)

})
