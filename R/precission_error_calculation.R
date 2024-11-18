#' Measurement Precision calculation
#'
#' Calculates the arithmetic differences between duplicate measurement pairs.
#' If the same measure have been measured more than twice, it computes the mean of the differences.
#' The differences are divided by sqrt(2) because they arise from imprecision in both
#' measurements, which adds quadratically, and the desired result is the precision
#' of one measurement.
#' @param tb A dataframe with id (unique ID value per duplicated sample), sensor (type of sensors) and measurements (the results of the measurements) columns
#' @return A tibble grouped by sensor that contains a dataframe with the next three columns 1.difference between duplicates, 2.mean values of each duplicate and 3.number of samples of the duplicates (commonly there will be just 2 sample as the name duplicate says, but if there are triplicates or more, the program will work as well is why it is important to know the n of samples) and RMS, MAD precision values and the number of samples that have been used to calculate these values (n_precision).
#' @examples
#' tb <- data.frame(id = c("m1", "m1", "m1", "m1", "m2", "m2", "m2", "m2"),
#'                  sensor = c("s1", "s1", "s2", "s2", "s1", "s1", "s2", "s2"),
#'                  measurement = c(10, 9, 2, 1, 15, 17, 4, 3))
#' tb_precision <- measure_precision(tb);
#' @export

measure_precision <- function(tb){

  defined <- ls()
  rlog::log_trace("Starting Measure Precission calculation")
  passed <- names(as.list(match.call())[-1])

  if (!defined %in% passed) {
    cli::cli_h1("Measurement's Precission calculation Outputs")
    cli::cli_li("{.strong Measurement difference Plot:} Differences between the sensor duplicate measurements are visualized in absolute and relative values (difference / mean)")
    cli::cli_li("{.strong Precission calculation:} (For patients with duplicate values, result = error in %) The relative differences of the duplicates are used to calculate the RMS and MAD precission for each of the sensors.")
    cli::cli_h1("Measurement's Precission calculation Inputs")
    cli::cli_li("{.strong Dataframe:} Dataframe with id (unique id to identify duplicates), sensor (the type of sensor used to measure) and measurement (measurement results, such as the concentration values of CO2 in the air) columns")
    return(invisible())
  }

  required_columns <- c("id", "sensor", "measurement")

  # Check if required columns are present in the dataframe
  missing_columns <- setdiff(required_columns, colnames(tb))

  # If there are missing columns, show a message
  if (length(missing_columns) > 0) {
    missing_msg <- paste("The following required columns are missing:", paste(missing_columns, collapse = ", "))
    stop(missing_msg) # Stop execution with an error message
  }

  if (sum(is.na(tb$id)) > 0) {
    missing_msg <- paste("Your dataframe have", sum(is.na(tb$id)), "NAs in id column")
    stop(missing_msg) # Stop execution with an error message
  }

  require(dplyr)
  require(purrr)
  require(tidyr)
  require(ggplot2)
  require(gridExtra)

  # FUNCTION STARTS HERE-------------------------------------------------------

  # Covert into factor the sensor column
  tb$sensor <- as.factor(tb$sensor)

  # Eliminate any NA value
  tb <- tb %>%
    dplyr::filter(!is.na(sensor))%>%
    dplyr::filter(!is.na(measurement))

  # Calculate the difference of each sensor measurement
  tb_diff <- tb %>%
    dplyr::group_by(id, sensor) %>%
    dplyr::summarise(
      diff = mean_pairwise_diff(measurement)[1],
      mn = mean(measurement),
      n = mean_pairwise_diff(measurement)[2],
    )

  # Make sure that the the length of tb_diff is correct
  # Make sure that there are not NA values in the dataframe
  # In the case there are NA values delete them
  # Report how many NA values and how many non-NA values there are

  # Plot the absolute and relative difference -------------------
  a <- ggplot(tb_diff, aes(x=mn, y=diff))+
    geom_point(aes(colour = sensor, shape = sensor))+
    ggtitle("Absolute difference")+
    labs(x = "Quantity mean", y = "Quantity diff")

  b <- ggplot(tb_diff, aes(x=mn, y=diff/mn))+
    geom_point(aes(colour = sensor, shape = sensor))+
    # ylim(c(0, 1))+
    ggtitle("Relative difference")+
    labs(x = "Quantity mean", y = "Relative Quantity diference")

  grid.arrange(a , b,  ncol=2, nrow=1)

  # PRECISSION CALCULOUS----------------------------------------

  # Create a tibble with the separate dataframes of the different miRNAs
  tb_diff2 <- tb_diff%>%
    dplyr::group_by(sensor)%>%
    tidyr::nest()%>%
    dplyr::ungroup()

  # Apply the MAD and RMS precision calculous
  tb_diff2 <- tb_diff2 %>%
    dplyr::mutate(RMS = purrr::map(data,RMS_precision))%>%
    dplyr::mutate(MAD = purrr::map(data,MAD_precision))%>%
    dplyr::mutate(n_precission = purrr::map(data, function(x) nrow(x)))

  return(tb_diff2)
}

#' Error Propagation calculation
#'
#' Calculates propagated errors for measurements based on an input formula, taking into account relative errors for each measurement component.
#' This function calculates the combined error for a new measurement by finding how each input's error affects the result.
#' It assumes that the errors in each input measurement are small and do not depend on each other.
#'
#' @param tb A tibble containing the measurements to be used in the error propagation calculation. The tibble should include
#' columns for `sensor` (the type of sensor), `err` (error of the sensor in %) and `data` (dataframe with measurement values by sensor with columns for `id` (sample unique id), `mn` (mean of the duplicated measurements) and `n` (number of identical sample, in case of duplicates 2 but sometimes some of them have more than 2 copies)).
#' This parameter is the output of the function `measure_precision` (it just need to change the name of the column of the error that wanted to be used to `err`).
#' @param fun An R expression representing the formula for calculating a new measurement based on input variables (sensor data).
#' This formula should specify the relationship between measurements in `tb` and the desired output measurement.
#' @return A dataframe with calculated measurement values (`measure`) and their propagated errors (`error`) for each measurement ID.
#' @examples
#' # Example dataframe
#' tb <- data.frame(sensor = c("s1", "s2"), data = list(data.frame(id = c("m1", "m2"), mn = c(10, 15), err = c(0.5, 0.6), n = c(2, 2))))
#' # Example formula
#' fun <- expression(s1 + s2)
#' # Calculate propagated errors
#' result <- error_propagation(tb, fun)
#' @export

error_propagation <- function (tb, fun){

  # save the variable names from the function
  vars <- all.vars(fun)
  # Sort the var names ascending alphabetically
  vars <- sort(vars, decreasing = FALSE)

  df_measure <- tb %>%
    dplyr::select(sensor, data)%>%
    # Remove diff column from the data dataframe
    dplyr::mutate(data = purrr::map(data, ~ .x %>% dplyr::select(id, mn))) %>%
    tidyr::unnest(data) %>%
    tidyr::pivot_wider(names_from = sensor, values_from = mn)%>%
    dplyr::arrange(id)

  # Add the variable nsensor to later save the information of the counts used for
  # the calculation of the mean value (This will be important in order to
  # calculate the real variability of that measurement)
  tb$nsensor <- paste0("n_", tb$sensor)


  # Save the n_miRNA value
  df_n <- tb%>%
    dplyr::select(nsensor, data)%>%
    # Remove diff and n column from the data dataframe
    dplyr::mutate(data = purrr::map(data, ~ .x %>% dplyr::select(id, n))) %>%
    tidyr::unnest(data) %>%
    tidyr::pivot_wider(names_from = nsensor, values_from = n)%>%
    dplyr::arrange(id)%>%
    dplyr::select(-id)

  # Merge both dataframes
  df_measure <- cbind(df_measure, df_n)

  derivatives <- lapply(vars, function(var) D(fun, var))

  # Ensure that the position of the sensors in the tibble is the same as in the "vars"
  tb <- tb %>%
    dplyr::arrange(sensor)

  # Get the relative error values for each var
  err <- as.numeric(tb$err[tb$sensor %in% vars])/(100)

  # create the Expression for error propagation
  error_prop <- 0

  for (i in seq_along(derivatives)) {
    # Create the expression part by part

    term <- substitute(der^2 * (error/sqrt(n_val)*val)^2, list(der = derivatives[[i]], error = err[i], val = str2lang(vars[i]), n_val = str2lang(paste0("n_",vars[i]))))
    error_prop <- substitute(a + b, list(a = error_prop, b = term))
  }

  error_prop <- substitute(sqrt(x), list(x = error_prop))

  apply_equation <- function(...) {

    # It is neccessary to have a df with the function variable names as columns

    # Capture the inputs in a list
    inputs <- list(...)

    # Evaluate the expression, mapping inputs to x1, x2, ...
    res <- eval(fun[[1]], envir = inputs)

    err <- eval(do.call(substitute, list(error_prop, inputs)))

    return(c(res, err))
  }


  # Aplly the function in the dataframe
  results <- do.call(mapply, c(FUN = apply_equation, as.list(df_measure)))
  rownames(results) <- c("measure", "error")

  # Result have 3 rows, one with the mean value and then the RMS and MAD values.
  # The column names are the IDs

  # Reshape the results
  results2 <- as.data.frame(results) %>%
    tibble::rownames_to_column("metric") %>%   # Move row names to a column called "metric"
    tidyr::pivot_longer(-metric, names_to = "ID", values_to = "value") %>%
    # Pivot to long format, the column names goes to ID, and values are distributed.
    # The measure type is defined in metric column (measure or error)
    # And the value in value column
    tidyr::pivot_wider(names_from = metric, values_from = value)

  return(results2)
}

