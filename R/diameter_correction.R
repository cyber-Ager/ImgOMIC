#' Correct diameter values for a single-patient dataset
#'
#' Adjusts measured artery diameters by accounting for measurement error
#' based on imaging modality (Ultrasound or CT). It generates multiple
#' candidate diameters with associated probabilities derived from a normal
#' distribution of measurement error, then uses a Viterbi-like algorithm to
#' find all most-probable non-decreasing diameter progression curves.
#'
#' @param df_patient A data.frame with columns:
#'   - `Date`: Numeric time points (e.g., decimal years)
#'   - `Diam`: Measured diameters
#'   - `CT`: Indicator (1=CT, 0=Ultrasound)
#' @param sdUS Numeric. Standard deviation of Ultrasound diameter measurements.
#'   Default: 2.
#' @param sdCT Numeric. Standard deviation of CT-scan diameter measurements.
#'   Default: 1.
#' @param sp Integer. Number of discretization points for the error distribution.
#'   Default: 4.
#' @param dlim_inf Numeric. Minimum allowable growth rate
#'   `(Diam_t - Diam_{t-1}) / (Date_t - Date_{t-1})`. Default: 0.
#' @param dlim_sup Numeric. Maximum allowable growth rate
#'   `(Diam_t - Diam_{t-1}) / (Date_t - Date_{t-1})`. Default: 50.
#'
#' @return A list containing:
#'   - `curves`: A list of the most-probable diameter progression vectors.
#'   - `max_prob`: The maximum joint probability among these curves.
#'   - `num_curves`: Number of curves sharing that maximum probability.
#'
#' @examples
#' data <- data.frame(
#'   ID = factor(c("1","1","1","1","1","1","1","1","1",
#'                 "2","2","2","2")),
#'   Date = c(2015, 2016, 2017, 2020, 2021, 2024, 2025, 2026, 2027,
#'            2013, 2016, 2018, 2019),
#'   Diam = c(30, 35, 34, 53, 50, 52, 55, 55, 57, 38, 42, 50, 53),
#'   CT   = factor(c(1,0,0,1,0,0,1,0,1,0,0,0,1))
#' )
#' df_patient <- subset(data, ID == "1")
#' res <- correct_diameter_single(df_patient)
#' @export
correct_diameter_single <- function(df_patient,
                                    sdUS = 3.5, sdCT = 1.9, sp = 4,
                                    dlim_sup = 1000, dlim_inf = 0) {

  # Check if input is empty
  if (nrow(df_patient) == 0) {
    stop("Input data frame is empty.")
  }

  # Check if input is empty
  if (nrow(df_patient) == 1) {
    stop("Input data-frame needs more than one value to be computed.")
  }

  if (any(is.na(df_patient$Date) |
          is.na(df_patient$CT) |
          is.na(df_patient$Diam))) {
    stop("Input data-frame has NA values in some of the key columns: Date, CT, Diam")
  }

  # Diam need to be numeric
  if (!is.numeric(df_patient$Diam)){

    stop("Diam column need to be numeric")

  }

  # CT need to have 0 and 1 values only
  if (!all(unique(as.character(df_patient$CT)) %in% c("0" , "1"))){

    stop("CT column can only handle 1 and 0 values")

  }

  if (length(unique(df_patient$Date))!= length(df_patient$Date)){

    stop("Date column can not have duplicated values")

  }

  # Arrange rows by the Date
  df_patient <- df_patient %>%
    dplyr::arrange(Date)

  # Check required columns
  required_cols <- c("Date", "Diam", "CT")
  missing_cols <- setdiff(required_cols, names(df_patient))
  if (length(missing_cols) > 0) {
    stop("Input data frame is missing required columns: ", paste(missing_cols, collapse = ", "))
  }

  # Calculate probability distribution limits
  limUS <- 3 * sdUS  # ±3 SD for Ultrasound
  limCT <- 3 * sdCT  # ±3 SD for CT

  # Discretize error distributions
  pUSx <- seq(-limUS, limUS, length.out = sp + 1)
  pCTx <- seq(-limCT, limCT, length.out = sp + 1)
  pUSy <- dnorm(pUSx, 0, sdUS); pUSy <- pUSy / sum(pUSy)
  pCTy <- dnorm(pCTx, 0, sdCT); pCTy <- pCTy / sum(pCTy)

  # Build fuzzy distributions for each time point
  dist_list <- purrr::map2(df_patient$Diam, df_patient$CT, function(diam, ct) {
    if (ct == 1) list(values = diam + pCTx, probs = pCTy)
    else         list(values = diam + pUSx, probs = pUSy)
  })

  diam_distr <- purrr::map(dist_list, "values")
  diam_prob  <- purrr::map(dist_list, "probs")
  n <- length(diam_distr)
  paths <- list()

  # Initialize at t=1
  paths[[1]] <- data.frame(
    value    = diam_distr[[1]],
    prob     = diam_prob[[1]],
    prev_idx = I(rep(list(NA), length(diam_distr[[1]])))
  )

  # Recursively build most probable paths
  for (t in 2:n) {
    prev_states <- paths[[t-1]]
    curr_vals   <- diam_distr[[t]]
    curr_probs  <- diam_prob[[t]]
    dt <- df_patient$Date[t] - df_patient$Date[t-1]

    value_list <- c(); prob_list <- c(); prev_list <- list()
    for (j in seq_along(curr_vals)) {
      val_j <- curr_vals[j]
      prob_j <- curr_probs[j]

      rate <- (val_j - prev_states$value) / dt
      valid_prev <- which(rate >= dlim_inf & rate <= dlim_sup)
      if (!length(valid_prev)) next

      acc_probs <- prev_states$prob[valid_prev] * prob_j
      max_p    <- max(acc_probs)
      best_prev<- valid_prev[which(acc_probs == max_p)]

      value_list <- c(value_list, val_j)
      prob_list  <- c(prob_list,  max_p)
      prev_list  <- c(prev_list,  list(best_prev))
    }

    if (!length(value_list)) {
      warning("No valid transitions at t = ", t)
      return(NULL)
    }

    paths[[t]] <- data.frame(
      value    = value_list,
      prob     = prob_list,
      prev_idx = I(prev_list)
    )
  }

  # Backtrack all most probable curves
  last_probs <- paths[[n]]$prob
  max_prob   <- max(last_probs)
  final_idxs <- which(last_probs == max_prob)

  build_paths <- function(t, idx) {
    if (t == 1) return(list(paths[[1]]$value[idx]))
    prev_idxs <- paths[[t]]$prev_idx[[idx]]
    unlist(lapply(prev_idxs, function(p) {
      lapply(build_paths(t-1, p), function(cur) c(cur, paths[[t]]$value[idx]))
    }), recursive = FALSE)
  }

  all_best_paths <- unlist(lapply(final_idxs, function(i) build_paths(n, i)), recursive = FALSE)

  list(curves     = all_best_paths,
       dates      = df_patient$Date,
       max_prob   = max_prob,
       num_curves = length(all_best_paths))
}

#' Correct diameter values across multiple patients
#'
#' Applies the single-patient correction to each patient in a dataset.
#' Converts Date to decimal years if needed, filters out rows with missing
#' Diam, CT, or Date, and replaces Diam with the most-probable curve (or the
#' average when multiple curves tie). Adds a logical `corrected` flag.
#'
#' @param data A data.frame with columns ID, Date, Diam, CT.
#' @param sdUS Numeric. SD for Ultrasound measures. Default: 2.
#' @param sdCT Numeric. SD for CT measures. Default: 1.
#' @param sp Integer. Number of discretization points. Default: 4.
#' @param dlim_inf Numeric. Minimum growth rate. Default: 0.
#' @param dlim_sup Numeric. Maximum growth rate. Default: 50.
#'
#' @return A data.frame with the same rows as `data`, updated `Diam`, and
#' a new `corrected` column indicating which values were adjusted.
#'
#' @examples
#' data <- data.frame(
#'   ID = factor(c("1","1","2","2")),
#'   Date = c(2015, 2016, 2015, 2017),
#'   Diam = c(30, 35, 40, 45),
#'   CT   = factor(c(1, 0, 1, 0))
#' )
#' res_all <- correct_diameters_all(data)
#' @export
correct_diameters_all <- function(data, sdUS = 3.5, sdCT = 1.9, sp = 4,
                                  dlim_inf = 0, dlim_sup = 50) {

  # Check if input is a data frame and not empty
  if (!is.data.frame(data) || nrow(data) == 0) {
    stop("`data` must be a non-empty data frame.")
  }

  # Check required columns
  required_cols <- c("ID", "Date", "Diam", "CT")
  missing_cols <- setdiff(required_cols, names(data))
  if (length(missing_cols) > 0) {
    stop("Input data frame is missing required columns: ", paste(missing_cols, collapse = ", "))
  }

  # Convert ID to character to forget about ghost levels
  data$ID <- as.character(data$ID)

  # Save previous date
  pre_date <- data$Date

  # Converts date to numeric format
  if (inherits(data$Date, "Date")) {
    data <- dplyr::mutate(data, Date = decimal_date(Date))
  }

  data <- dplyr::mutate(data, CT = as.numeric(as.character(CT)))


  # map_dfr make any process in teh data and afterwards it bind the rows and cols
  out <- data %>%
    split(.$ID) %>%
    purrr::map_dfr(function(df_patient) {

      # Add column
      df_patient$corrected <- FALSE

      # Check valid ids
      valid_idx <- which(!is.na(df_patient$Diam) &
                           !is.na(df_patient$CT)   &
                           !is.na(df_patient$Date))

      # Get valid data frame
      df_valid  <- df_patient[valid_idx, ]

      # Process it if it has good length
      if (nrow(df_valid) >= 2) {
        res <- correct_diameter_single(df_valid, sdUS, sdCT, sp,
                                       dlim_inf = dlim_inf,
                                       dlim_sup = dlim_sup)
        if (!is.null(res$curves)) {
          if (res$num_curves == 1) {
            best_curve <- res$curves[[1]]
          } else {
            best_curve <- colMeans(do.call(rbind, res$curves))
          }
          # Modify the previous dataframe
          # Taking into account the previous order of the df
          df_patient$Diam[valid_idx]      <- best_curve[order(df_patient$Date[valid_idx])]
          df_patient$corrected[valid_idx] <- TRUE
        }
      }
      df_patient
    })
  # Save previous date format
  out$Date <- pre_date
  out
}
