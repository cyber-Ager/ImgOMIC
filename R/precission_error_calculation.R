#' Measurement Precision calculation
#'
#' Calculates the arithmetic differences between duplicate measurement pairs.
#' If the same measure have been measured more than twice, it computes the mean of the differences.
#' The differences are divided by sqrt(2) because they arise from imprecision in both
#' measurement, which adds quadratically, and the desired result is the precision
#' of one measurement.
#' @param tb A dataframe with `ID` (unique ID value per duplicated sample), `sensor` (type of sensors) and `measurements` (the results of the measurements) columns
#' @param img (LOGICAL element) It determines if you wants plots or not. It is TRUE as default.
#' @return A tibble grouped by sensor that contains a dataframe with the next three columns 1.difference between duplicates, 2.mean values of each duplicate and 3.number of samples of the duplicates (commonly there will be just 2 sample as the name duplicate says, but if there are triplicates or more, the program will work as well is why it is important to know the n of samples) and RMS, MAD precision values and the number of samples that have been used to calculate these values (n_precision).
#' @examples
#' tb <- data.frame(ID = c("m1", "m1", "m1", "m1", "m2", "m2", "m2", "m2"),
#'                  sensor = c("s1", "s1", "s2", "s2", "s1", "s1", "s2", "s2"),
#'                  measurement = c(10, 9, 2, 1, 15, 17, 4, 3))
#' tb_precision <- measure_precision(tb);
#' @export

measure_precision <- function(tb, img = TRUE){

  required_columns <- c("ID", "sensor", "measurement")

  # Check if required columns are present in the dataframe
  missing_columns <- setdiff(required_columns, colnames(tb))

  # If there are missing columns, show a message
  if (length(missing_columns) > 0) {
    missing_msg <- paste("The following required columns are missing:", paste(missing_columns, collapse = ", "))
    stop(missing_msg) # Stop execution with an error message
  }

  if (sum(is.na(tb$ID)) > 0) {
    missing_msg <- paste("Your dataframe have", sum(is.na(tb$ID)), "NAs in ID column")
    stop(missing_msg) # Stop execution with an error message
  }

  # Install the required packages
  list.of.packages <- c("ggplot2", "dplyr", "purrr", "tidyr", "gridExtra")
  new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
  if(length(new.packages)) install.packages(new.packages)

  library(dplyr)

  # FUNCTION STARTS HERE-------------------------------------------------------

  # Covert into factor the sensor column
  tb$sensor <- as.factor(tb$sensor)

  # Ensure that the column measurement is numeric and eliminate anything that is not a number
  tb$measurement <- as.numeric(tb$measurement)

  # Eliminate any NA value
  tb <- tb %>%
    dplyr::filter(!is.na(sensor))%>%
    dplyr::filter(!is.na(measurement))

  if(length(tb$ID) == 0) {

    stop("All IDs have NA values in measurement OR sensor columns")

  }

  # Calculate the difference of each sensor measurement
  tb_diff <- tb %>%
    dplyr::group_by(ID, sensor) %>%
    dplyr::summarise(
      diff = mean_pairwise_diff(measurement)[1],
      mn = mean(measurement),
      n = mean_pairwise_diff(measurement)[2],
    )

  if (all(is.na(tb_diff$diff))) {

    stop("There are no duplicates in your dataframe")

  }

  # Count the single values

  non_duplicates <- sum(tb_diff$n == 1)

  if (non_duplicates > 0){
    message(paste(non_duplicates, "measurements have been removed since they do not have duplicate values"))
  }
  # Eliminate the NO-duplicate info
  tb_diff <- tb_diff %>%
    dplyr::filter(!is.na(diff))

  # Plot the absolute and relative difference -------------------
  if (img) {
    library(ggplot2)

    a <- ggplot(tb_diff, aes(x=mn, y=diff))+
      geom_point(aes(colour = sensor, shape = sensor))+
      ggtitle("Absolute difference")+
      labs(x = "Quantity mean", y = "Quantity diff")

    b <- ggplot(tb_diff, aes(x=mn, y=diff/mn))+
      geom_point(aes(colour = sensor, shape = sensor))+
      # ylim(c(0, 1))+
      ggtitle("Relative difference")+
      labs(x = "Quantity mean", y = "Relative Quantity diference")

    gridExtra::grid.arrange(a , b,  ncol=2, nrow=1)
  }

  # PRECISION CALCULOUS----------------------------------------

  # Create a tibble with the separate dataframes of the different miRNAs
  tb_diff2 <- tb_diff%>%
    dplyr::group_by(sensor)%>%
    tidyr::nest()%>%
    dplyr::ungroup()

  # Apply the MAD and RMS precision calculous
  tb_diff2 <- tb_diff2 %>%
    dplyr::mutate(RMS = purrr::map(data,RMS_precision))%>%
    dplyr::mutate(MAD = purrr::map(data,MAD_precision))%>%
    dplyr::mutate(n_precision = purrr::map(data, function(x) nrow(x)))

  return(tb_diff2)
}

#' Error Propagation calculation
#'
#' Calculates propagated errors for measurements based on an input formula, taking into account relative errors for each measurement component.
#' This function calculates the combined error for a new measurement by finding how each input's error affects the result.
#' It assumes that the errors in each input measurement are small and do not depend on each other.
#'
#' @param tb A tibble containing the measurements to be used in the error propagation calculation. The tibble should include
#' columns for `sensor` (the type of sensor), `err` (error of the sensor in %) and `data` (dataframe with measurement values by sensor with columns for `ID` (sample unique ID), `mn` (mean of the duplicated measurements) and `n` (number of identical sample, in case of duplicates 2 but sometimes some of them have more than 2 copies)).
#' This parameter is the output of the function `measure_precision` (it just need to change the name of the column of the error that wanted to be used to `err`).
#' @param fun An R expression representing the formula for calculating a new measurement based on input variables (sensor data).
#' This formula should specify the relationship between measurements in `tb` and the desired output measurement.
#' @return A dataframe with calculated measurement values (`measure`) and their propagated errors (`error`) for each measurement ID.
#' @examples
#' # Example dataframe
#' tb <- data.frame(sensor = c("s1", "s2"), data = list(data.frame(ID = c("m1", "m2"), mn = c(10, 15), err = c(0.5, 0.6), n = c(2, 2))))
#' # Example formula
#' fun <- expression(s1 + s2)
#' # Calculate propagated errors
#' result <- error_propagation(tb, fun)
#' @export

error_propagation <- function (tb, fun = NULL){

  if(is.null(fun)){

    stop("fun is missing. It is neccessary to provide the fun expression.")

  }

  if(class(fun) != "expression"){

    stop("fun need to be expression type")

  }

  required_columns <- c("sensor", "err", "data")

  # Check if required columns are present in the dataframe
  missing_columns <- setdiff(required_columns, colnames(tb))

  # If there are missing columns, show a message
  if (length(missing_columns) > 0) {
    missing_msg <- paste("The following required columns are missing:", paste(missing_columns, collapse = ", "))
    stop(missing_msg) # Stop execution with an error message
  }

  required_columns <- c("ID", "mn", "n")

  # Check if required columns are present in the dataframe
  missing_columns <- setdiff(required_columns, colnames(tb$data[[1]]))

  # If there are missing columns in the data dataframe, show a message
  if (length(missing_columns) > 0) {
    missing_msg <- paste("The following required columns are missing:", paste(missing_columns, collapse = ", "))
    stop(missing_msg) # Stop execution with an error message
  }

  # Install the required packages
  list.of.packages <- c("dplyr", "purrr", "tidyr", "tibble")
  new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
  if(length(new.packages)) install.packages(new.packages)

  library(dplyr)
  # FUNCTION STARTS HERE-----------------------------------------------------

  # save the variable names from the function
  vars <- all.vars(fun)
  # Sort the var names ascending alphabetically
  vars <- sort(vars, decreasing = FALSE)

  # Check if there is some variable in the function that is not defined in the tb
  missing_vars <- setdiff(vars, tb$sensor)
  if (length(missing_vars) > 0) {
    missing_msg <- paste("The following variables defined in the equation are not shown in the input tibble:", paste(missing_vars, collapse = ", "))
    stop(missing_msg) # Stop execution with an error message
  }

  error_check <- as.double(tb$err[tb$sensor == rep(vars, length.out = length(tb$sensor))])

  # Check if there is some error value that have not being defined
  if (sum(is.na(error_check)) > 0){

  sens_errors <- as.character(tb$sensor[is.na(error_check)])

    stop(paste("The following sensors have NA values in err column:", paste(sens_errors, collapse = ", ")))

  }

  # Ensure that mn and n are numeric values
  tb <- tb %>%
    dplyr::mutate(data = purrr::map(data, ~ .x %>% dplyr::mutate(
      mn = as.numeric(mn),
      n = as.numeric(n)
    )))


  # Simplify the tibble to save the mn value of each sensor divided in columns and ordered in rows by ID.
  df_measure <- tb %>%
    dplyr::select(sensor, data)%>%
    dplyr::mutate(data = purrr::map(data, ~ .x %>% dplyr::select(ID, mn))) %>%
    tidyr::unnest(data) %>%
    tidyr::pivot_wider(names_from = sensor, values_from = mn)%>%
    dplyr::arrange(ID)

  # Delete the rows where there was an NA value and inform about it
  df_measure2 <- df_measure %>%
    dplyr::filter(if_all(everything(), ~ !is.na(.)))

  df_dif <- nrow(df_measure) - nrow(df_measure2)

  if (df_dif > 0) {

    message(paste(df_dif, "samples have been excluded from the calulation due to missing values in sensor data."))

    df_measure <- df_measure2

  }

  if(nrow(df_measure) == 0){

    stop("0 samples with ALL the needed measurements. Impossible to calculate the new measurement in any sample.")

  }

  # Add the variable nsensor to later save the information of the counts used for
  # the calculation of the mean value (This will be important in order to
  # calculate the real variability of that measurement)
  tb$nsensor <- paste0("n_", tb$sensor)


  # Save the n_sensor value in the same way as the mn was saved before
  df_n <- tb%>%
    dplyr::select(nsensor, data)%>%
    dplyr::mutate(data = purrr::map(data, ~ .x %>% dplyr::select(ID, n))) %>%
    tidyr::unnest(data) %>%
    tidyr::pivot_wider(names_from = nsensor, values_from = n)%>%
    dplyr::arrange(ID)%>%
    dplyr::filter(if_all(everything(), ~ !is.na(.)))%>%
    dplyr::select(-ID)


  # Merge both dataframes
  df_measure <- cbind(df_measure, df_n)

  # Define the partial derivative of fun respect each sensor
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

  # Result have 2 rows containing the result of the function and the propagated error
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

