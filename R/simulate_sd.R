#' Monte Carlo simulation of patient diameter statistic under measurement uncertainty
#'
#' Performs a Monte Carlo simulation to estimate the mean and standard deviation
#' of a user-defined statistic (by default, the slope of a linear regression of
#' diameter vs. time), accounting for measurement error specific to imaging modality
#' (Ultrasound or CT). Each simulation perturbs observed diameters according to
#' assumed normal measurement error distributions and evaluates the chosen statistic.
#'
#' @param df_patient A `data.frame` with columns:
#'   - `Date`: Numeric time points (e.g., years or days)
#'   - `Diam`: Measured diameters (numeric)
#'   - `CT`: Indicator of imaging modality (1 = CT, 0 = Ultrasound; may be factor or numeric)
#' @param FUN A function of the form `FUN(Date, Diam)` returning a numeric scalar.
#'   Defines the statistic to be computed for each simulated dataset.
#'   If `NULL`, defaults to the slope from a simple linear regression of
#'   `Diam ~ Date`.
#' @param sdUS Numeric. Standard deviation of Ultrasound measurement error. Notice that you need to be consistent with the measurement units with Diam.
#'   Default: 3.5.
#' @param sdCT Numeric. Standard deviation of CT measurement error. Notice that you need to be consistent with the measurement units with Diam.
#'   Default: 1.9.
#' @param n_sim Integer. Number of Monte Carlo simulation iterations.
#'   Default: 100.
#' @param seed Optional integer. Random seed for reproducibility.
#'   Default: `NULL`.
#' @param min_diam Numeric. Minimum physically plausible diameter.
#'   Values below this are truncated. Default: 15.
#' @param max_diam Numeric. Maximum physically plausible diameter.
#'   Values above this are truncated. Default: 150.
#'
#' @return A list with elements:
#'   - `mean`: Mean of the simulated statistic across valid simulations.
#'   - `sd`: Standard deviation of the simulated statistic.
#'   - `n_valid`: Number of valid simulations.
#'
#' @details
#' For each simulation, diameters are perturbed using a normal distribution with
#' standard deviation determined by the imaging modality (`sdUS` or `sdCT`).
#' The provided function `FUN` is evaluated on each perturbed dataset. Simulations
#' producing non-finite results are excluded from the final summary.
#'
#' @examples
#' df <- data.frame(
#'   Date = c(2015, 2016, 2017, 2018),
#'   Diam = c(30, 34, 33, 36),
#'   CT   = c(0, 1, 0, 1)
#' )
#'
#' # Default usage: slope of linear regression of Diam ~ Date
#' res1 <- simulate_function_sd(df)
#'
#' # Custom statistic: mean diameter
#' res2 <- simulate_function_sd(df, FUN = function(d, x) mean(x))
#'
#' @export
simulate_function_sd <- function(df_patient,
                                            FUN = NULL,
                                            sdUS = 3.5,
                                            sdCT = 1.9,
                                            n_sim = 100,
                                            seed = NULL,
                                            min_diam = 15,
                                            max_diam = 150) {
  # df_patient: data.frame with Date (numeric), Diam (numeric), CT (0/1 or factor)
  # FUN: function(Date, Diam) -> numeric scalar. If NULL, defaults to slope of linear regression.
  # sdUS, sdCT: measurement standard deviations for US and CT respectively.
  # n_sim: number of Monte Carlo simulations
  # seed: optional integer for reproducibility
  # keep_samples: if TRUE returns all simulated statistic values (may be large)
  # min_diam: kimits for negative diameters
  # max_diam: limits a maximum diameter

  if (!is.data.frame(df_patient)) stop("df_patient must be a data.frame.")
  required_cols <- c("Date","Diam","CT")
  if (!all(required_cols %in% names(df_patient))) {
    stop("df_patient must contain columns: Date, Diam, CT")
  }
  # if (nrow(df_patient) < 2) stop("Need at least two time points.")
  # if (any(is.na(df_patient$Date) | is.na(df_patient$Diam) | is.na(df_patient$CT))) {
  #   stop("NA values found in Date, Diam or CT columns.")
  # }
  if (!is.numeric(df_patient$Date) || !is.numeric(df_patient$Diam)) {
    stop("Date and Diam must be numeric.")
  }
  # Normalize CT to numeric 0/1
  df_patient$CT <- as.numeric(as.character(df_patient$CT))
  if (!all(df_patient$CT %in% c(0,1,NA))) stop("CT column must contain only 0 and 1 (or factors convertible to that).")

  # default FUN: slope of linear regression (can be negative, large, etc.)
  if (is.null(FUN)) {
    FUN <- function(dates, diams) {
      # return slope (coef of date). If fit fails (e.g. constant diams) return NA
      if (length(unique(diams)) <= 1 || length(unique(dates)) <= 1) return(NA_real_)
      fit <- try(stats::lm(diams ~ dates), silent = TRUE)
      if (inherits(fit, "try-error")) return(NA_real_)
      coef(fit)[["dates"]]
    }
  }
  else {
    if (!is.function(FUN)) stop("FUN must be a function.")
    if (length(formals(FUN)) != 2) {
      stop("FUN must take exactly two arguments: FUN(dates, diams).")
    }
  }




  if (!is.null(seed)) set.seed(seed)

  # Number of entries
  n <- nrow(df_patient)
  dates <- df_patient$Date
  measured <- df_patient$Diam
  ctflags <- as.integer(df_patient$CT)

  stats_vec <- numeric(n_sim)

  # Vector with SD of each diameter measurement
  sds <- ifelse(ctflags == 1, sdCT, sdUS)

  for (i in seq_len(n_sim)) {
    # sample plausible true diameters for this iteration
    sampled <- rnorm(n, mean = measured, sd = sds)
    # Ensure that the values are in physically possible range
    sampled <- pmax(sampled, min_diam)
    sampled <- pmin(sampled, max_diam)

    # evaluate user function robustly
    val <- try(FUN(dates, sampled), silent = TRUE)
    if (inherits(val, "try-error") || length(val) != 1 || !is.finite(val)) {
      stats_vec[i] <- NA_real_
    } else {
      stats_vec[i] <- as.numeric(val)
    }
  }

  # remove NA entries (optional: warn if many NA)
  na_count <- sum(is.na(stats_vec))

  # In case there is some NAs print the next error
  if (na_count > 0) {
    warning(sprintf("%d/%d simulations produced NA for the statistic (they will be ignored).", na_count, n_sim))
  }

  stats_valid <- stats_vec[is.finite(stats_vec)]
  n_valid <- length(stats_valid)
  if (n_valid == 0) stop("All simulations failed to produce a finite statistic.")

  mean_stat <- mean(stats_valid)
  sd_stat   <- stats::sd(stats_valid)


  res <- list(
    mean = mean_stat,
    sd = sd_stat,
    n_valid = n_valid
  )
  return(res)
}
