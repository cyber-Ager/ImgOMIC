#' Correct diameter values of an artery from a single patient dataset.
#'
#' Takes into account the error of the diameter measures from an artery depending the technique that have been used for measuring it.
#' Generates additional diameter measures with an associated probability of being real based on the standard deviation of the measure.
#' Finds the most probable diameter progression curve using Viterbi algorithm.
#'
#' @param df_patient A data.frame with columns Date (numeric), Diam, CT.
#' @param sdUS Standard deviation of the diameter measured by Ultrasound, in the same units as the dataframe Diam values. 2 as default.
#' @param sdCT Standard deviation of the diameter measured by CT-scan in the same units as the dataframe Diam values. 1 as default.
#' @param sp Number of points to discretize the probability distribution. 4 as default.
#' @param dlim_inf Minimum allowed growth (Diam2-Diam1)/(Date2-Date1). 0 as default.
#' @param dlim_sup Maximum allowed growth (Diam2-Diam1)/(Date2-Date1). 50 as default.
#' @return A list with 3 elements. `curves` vector that contains the curves with the maximum probability of existence that fits the input conditions. `max_prob` The product of the probability of each point of the curve. `num_curves` number of curves with the maximum probability.
#' @examples
#'data <- data.frame(
#'  ID = as.factor(c("1", "1", "1", "1", "1", "1", "1", "1", "1", "2", "2", "2", "2")),
#'  Date = c(2015, 2016, 2017, 2020, 2021, 2024, 2025, 2026, 2027, 2013, 2016, 2018, 2019),
#'  Diam = c(30, 35, 34, 53, 50, 52, 55, 55, 57, 38, 42, 50, 53),
#'  CT = as.factor(c(1, 0, 0, 1, 0, 0, 1, 0, 1, 0, 0, 0, 1))
#')
#'
#' # Filter patients
#' df_patient <- data[ID==1,]
#' res <- correct_diameter_single(df_patient)
#' @export

# Get the corrected max curve
correct_diameter_single <- function(df_patient, sdUS = 2, sdCT = 1, sp =4,
                             dlim_sup = 1000, dlim_inf = 0) {


  limUS <- 3 * sdUS # Limits of the normal distribution in US (3xSD)
  limCT <- 3 * sdCT # Limits of the normal distribution in CT (3xSD)

  # x points of the normal distribution in US
  pUSx <- seq(-limUS, limUS, length.out = sp + 1)
  # x points of the normal distribution in CT
  pCTx <- seq(-limCT, limCT, length.out = sp + 1)
  # y points of the normal distribution in US ;  Ensure the sum is 1
  pUSy <- dnorm(pUSx, 0, sdUS); pUSy <- pUSy / sum(pUSy)
  # x points of the normal distribution in CT ;  Ensure the sum is 1
  pCTy <- dnorm(pCTx, 0, sdCT); pCTy <- pCTy / sum(pCTy)
  # Derivative limits to consider a curve real or not

  # Construcción de distribuciones borrosas
  dist_list <- purrr::map2(df_patient$Diam, df_patient$CT, function(diam, ct) {
    if (ct == 1) list(values = diam + pCTx, probs = pCTy)
    else         list(values = diam + pUSx, probs = pUSy)
  })

  diam_distr <- purrr::map(dist_list, "values")
  diam_prob  <- purrr::map(dist_list, "probs")
  n <- length(diam_distr)

  # paths: lista con data.frames para cada tiempo
  paths <- list()

  # Paso 1: inicialización
  paths[[1]] <- data.frame(
    value = diam_distr[[1]],
    prob = diam_prob[[1]],
    prev_idx = I(rep(list(NA), length(diam_distr[[1]])))  # lista de previos (en t = 1 → NA)
  )

  # Paso 2: recorrer el resto del tiempo
  for (t in 2:n) {
    prev_states <- paths[[t - 1]]
    curr_vals   <- diam_distr[[t]]
    curr_probs  <- diam_prob[[t]]

    # Preparar lista para el paso t
    value_list <- c()
    prob_list <- c()
    prev_list <- list()

    # Calculate time difference
    dt <- df_patient$Date[t] - df_patient$Date[t-1]

    for (j in seq_along(curr_vals)) {
      val_j <- curr_vals[j]
      prob_j <- curr_probs[j]



      # Compute rate of change (mm per year, for example)
      growth_rate <- (val_j - prev_states$value) / dt

      # Filter valid growth rates
      valid_prev <- which(growth_rate >= dlim_inf & growth_rate <= dlim_sup)

      if (length(valid_prev) == 0) next

      acc_probs <- prev_states$prob[valid_prev] * prob_j
      max_prob <- max(acc_probs)
      best_prev <- valid_prev[which(acc_probs == max_prob)]

      # Guardar estado actual
      value_list <- c(value_list, val_j)
      prob_list  <- c(prob_list, max_prob)
      prev_list  <- c(prev_list, list(best_prev))
    }

    if (length(value_list) == 0) {
      warning("No valid transitions at t = ", t)
      return(NULL)
    }

    paths[[t]] <- data.frame(
      value = value_list,
      prob = prob_list,
      prev_idx = I(prev_list)  # almacenar lista de índices anteriores válidos
    )
  }

  # Paso 3: reconstruir TODAS las curvas más probables
  last_probs <- paths[[n]]$prob
  max_prob <- max(last_probs)
  final_idxs <- which(last_probs == max_prob)

  # Backtracking recursivo
  build_paths <- function(t, idx) {
    if (t == 1) {
      return(list(paths[[1]]$value[idx]))
    }
    prev_idxs <- paths[[t]]$prev_idx[[idx]]
    sub_paths <- unlist(lapply(prev_idxs, function(pidx) {
      lapply(build_paths(t - 1, pidx), function(p) c(p, paths[[t]]$value[idx]))
    }), recursive = FALSE)
    return(sub_paths)
  }

  all_best_paths <- unlist(lapply(final_idxs, function(idx) {
    build_paths(n, idx)
  }), recursive = FALSE)

  return(list(
    curves = all_best_paths,
    max_prob = max_prob,
    num_curves = length(all_best_paths)
  ))
}


#' Correct diameter values of an artery from a multi-patient dataset.
#'
#' Takes into account the error of the diameter measures from an artery depending the technique that have been used for measuring it.
#' Generates additional diameter measures with an associated probability of being real based on the standard deviation of the measure.
#' Finds the most probable diameter progression curve using Viterbi algorithm and returns a dataframe with the corrected values.
#' The output dataset has an additional column where it specifies if the diameter progression curve have been corrected or not.
#'
#' @param df_patient A data.frame with columns ID, Date, Diam, CT.
#' @param sdUS Standard deviation of the diameter measured by Ultrasound, in the same units as the dataframe Diam values. 2 as default.
#' @param sdCT Standard deviation of the diameter measured by CT-scan in the same units as the dataframe Diam values. 1 as default.
#' @param sp Number of points to discretize the probability distribution. 4 as default.
#' @param dlim_inf Minimum allowed growth (Diam2-Diam1)/(Date2-Date1). 0 as default.
#' @param dlim_sup Maximum allowed growth (Diam2-Diam1)/(Date2-Date1). 50 as default.
#' @return A dataframe with the same number of rows as the input dataframe with corrected diameter values and a additional column that specifies if the measure have been corrected or not.
#' @examples
#'data <- data.frame(
#'  ID = as.factor(c("1", "1", "1", "1", "1", "1", "1", "1", "1", "2", "2", "2", "2")),
#'  Date = c(2015, 2016, 2017, 2020, 2021, 2024, 2025, 2026, 2027, 2013, 2016, 2018, 2019),
#'  Diam = c(30, 35, 34, 53, 50, 52, 55, 55, 57, 38, 42, 50, 53),
#'  CT = as.factor(c(1, 0, 0, 1, 0, 0, 1, 0, 1, 0, 0, 0, 1))
#')
#'
#' res <- correct_diameter_all(df_patient)
#' @export


correct_diameters_all <- function(data ,sdUS = 2, sdCT = 1, sp = 4,
                                  dlim_inf = 0, dlim_sup = 50) {

  # 1) Convert Date to decimal years if needed
  if (inherits(data$Date, "Date")) {
    data <- data %>% dplyr::mutate(Date = decimal_date(Date))
  }

  # 2) Ensure CT is numeric 0/1
  data <- data %>% dplyr::mutate(CT = as.numeric(as.character(CT)))

  # 3) Process per patient
  out <- data %>%
    split(.$ID) %>%
    purrr::map_dfr(function(df_patient) {

      # Prepare a flag column
      df_patient$corrected <- FALSE

      # Subset to non-NA Diam & CT
      valid_idx <- which(!is.na(df_patient$Diam) &
                           !is.na(df_patient$CT) &
                           !is.na(df_patient$Date))
      df_valid  <- df_patient[valid_idx, ]

      # Only try if at least 2 valid, and ≥3 total rows
      if (nrow(df_valid) >= 2) {
        res <- corrected_curves(df_valid, sdUS,sdCT,
                                dlim_inf = dlim_inf,
                                dlim_sup = dlim_sup)
        # If the function returned a result list with curves
        if (!is.null(res) && !is.null(res$curves)) {
          # 4) Pick or average
          if (res$num_curves == 1) {
            best_curve <- res$curves[[1]]
          } else {
            # build matrix (each row is one curve) and take mean per column
            mat <- do.call(rbind, res$curves)
            best_curve <- colMeans(mat)
          }
          # 5) Replace Diam and flag
          df_patient$Diam[valid_idx]      <- best_curve
          df_patient$corrected[valid_idx] <- TRUE
        }
      }

      return(df_patient)
    })

  return(out)
}
