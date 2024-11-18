
#' Mean pairwise calculation
#'
#' Calculates the arithmetic differences between duplicate measurement pairs.
#' If the same measure have been measured more than twice, it computes the mean of the differences.
#' The differences are divided by sqrt(2) because they arise from imprecision in both
#' measurements, which adds quadratically, and the desired result is the precision
#' of one measurement.
#' @param values The measurement duplicated values of the same sample.
#' @return The difference of the measurements divided by sqrt(2) and the number of samples used for the calculus.
#' @examples
#' diff1 <- mean_pairwise_diff(c(22, 26));
#' diff2 <- mean_pairwise_diff( c(10, 12, 15) );
#' @export

# Define a function to calculate the mean pairwise differences
mean_pairwise_diff <- function(values) {

  if (sum(is.na(values)) > 0) {

    stop("There are NA values in the mean_pairwise_diff calculous")

    }

  if (length(values) == 2) {
    return(c((values[1] - values[2])/sqrt(2), length(values)))
  }
  if (length(values) == 1){
    return(c(NA, length(values)))
  }
  else {
    # Calculate all pairwise combinations of differences
    combn_diffs <- combn(values, 2, function(x) abs(x[1] - x[2])/sqrt(2))
    return(c(mean(combn_diffs), length(values)))
  }
}

#' Root Mean Square precision calculation
#'
#' Calculates the root mean square precision of a dataframe that contains the
#' mean and difference values of measured duplicates.
#' @param df Dataframe containing the mean (mn) and difference (diff) columns
#' @return The Root Mean Square Precision value in %.
#' @examples
#' RMS1 <- RMS_precision(data.frame(mn = c(2 , 5 , 6), diff(0.1, 0.5, 0.2)));
#' @export

RMS_precision <- function(df){

  required_columns <- c("diff", "mn")

  # Check if required columns are present in the dataframe
  missing_columns <- setdiff(required_columns, colnames(df))

  # If there are missing columns, show a message
  if (length(missing_columns) > 0) {
    missing_msg <- paste("The following required columns are missing:", paste(missing_columns, collapse = ", "))
    stop(missing_msg) # Stop execution with an error message
  }

  # dif is the difference between samples and mean is the mean value

  D1 <- df$diff/df$mn
  n <- length(df$diff)

  RMS <- sqrt(sum(D1^2)/n)*100
  return(RMS)

}

#' Mean Absolute Difference (MAD) precision calculation
#'
#' Calculates the mean absolute difference precision of a dataframe that contains the
#' mean and difference values of measured duplicates.
#' @param df Dataframe containing the mean (mn) and difference (diff) columns
#' @return The MAD Precision value in %.
#' @examples
#' MAD1 <- MAD_precision(data.frame(mn = c(2 , 5 , 6), diff(0.1, 0.5, 0.2)));
#' @export

MAD_precision <- function(df){

  required_columns <- c("diff", "mn")

  # Check if required columns are present in the dataframe
  missing_columns <- setdiff(required_columns, colnames(df))

  # If there are missing columns, show a message
  if (length(missing_columns) > 0) {
    missing_msg <- paste("The following required columns are missing:", paste(missing_columns, collapse = ", "))
    stop(missing_msg) # Stop execution with an error message
  }

  # dif is the difference between samples and mean is the mean value

  D1 <- df$diff/df$mn
  n <- length(df$diff)

  MAD <- sqrt(pi/2)*sum(abs(D1))/n*100

  return(MAD)
}
