
# FUNCTION TO COMPUTE BOOTSTRAP CONFIDENCE INTERVALS --------------------------

#' Computes Bootstrap confidence intervals
#'
#' @param b An object of class 'boot', containing the output
#' of a bootstrap calculation.
#' @param conf A scalar or vector containing the confidence
#' level(s) of the required interval(s).
#' @param type A vector of character strings representing
#' the type of intervals required. The value should be any subset
#' of the values c('norm','basic', 'perc', 'bca').
#' @importFrom boot 'boot.ci'
#' @importFrom stats 'sd'
#'
#' @return An object of class boot
bootstats <- function(b, conf = conf, type = type) {
  p <- length(b$t0)
  lab <- c("original", "bias", "std.error", "low.ci", "high.ci")
  out <- as.data.frame(matrix(nrow = p, ncol = length(lab), dimnames = list(NULL, lab)))
  for (i in 1:p) {
    # original estimation, bias, standard deviation
    out[i, "original"] <- b$t0[i]
    out[i, "bias"] <- mean(b$t[, i]) - b$t0[i]
    out[i, "std.error"] <- sd(b$t[, i])
    # confidence interval
    if (type == "norm") {
      ci <- boot.ci(b, index = i, type = "norm", conf = conf)
      if (!is.null(ci)) {
        out[i, "low.ci"] <- ci$norm[2]
        out[i, "high.ci"] <- ci$norm[3]
      }
    } else if (type == "basic") {
      ci <- boot.ci(b, index = i, type = "basic", conf = conf)
      if (!is.null(ci)) {
        out[i, "low.ci"] <- ci$basic[4]
        out[i, "high.ci"] <- ci$basic[5]
      }
    } else if (type == "percent") {
      ci <- boot.ci(b, index = i, type = "perc", conf = conf)
      if (!is.null(ci)) {
        out[i, "low.ci"] <- ci$percent[4]
        out[i, "high.ci"] <- ci$percent[5]
      }
    } else if (type == "bca") {
      ci <- boot.ci(b, index = i, conf = conf)
      if (!is.null(ci)) {
        out[i, "low.ci"] <- ci$bca[4]
        out[i, "high.ci"] <- ci$bca[5]
      }
    }
  }
  return(out)
}

#' Computes bootstrap confidence intervals
#'
#' @param b A boot object with the computed Sobol' indices.
#' @param params A vector with the name of the model inputs.
#' @param type A vector of character strings representing
#' the type of intervals required. The value should be any subset
#' of the values c('norm','basic', 'perc', 'bca'). For more information,
#' check the function boot::boot.ci.
#' @param conf A scalar or vector containing the confidence
#' level(s) of the required interval(s).
#' @param second Boolean. If second == TRUE, it computes the confidence
#' intervals for second-order indices. Default is second == FALSE.
#' @param third Boolean. If third == TRUE, it computes the confidence
#' intervals for second-order indices. Default is third == FALSE.
#'
#' @return A data table.
sobol_ci_temp <- function(b, params, type, conf, second = FALSE, third = FALSE) {
  V1 <- NULL
  # Extract first and total order effects
  ext1 <- b[, .SD[1:length(params)]][, V1]
  Si.STi <- lapply(ext1, function(x) bootstats(x,
                                               type = type,
                                               conf = conf)) %>%
    data.table::rbindlist()
  if(second == TRUE) {
    # Extract second order effects
     ext2 <- b[, .SD[(length(params) + 1) :(length(params) +
                                              factorial(length(params)) /
                                              (factorial(2) * factorial(length(params) - 2)))]][, V1]
     Sij <- lapply(ext2, function(x) bootstats(x,
                                               type = type,
                                               conf = conf)) %>%
       data.table::rbindlist()
  } else {
    Sij <- NULL
  }
  if(third == TRUE) {
    ext3 <- b[, utils::tail(.SD, factorial(length(params)) /
                              (factorial(3) * factorial(length(params) - 3)))][, V1]
    Sijk <- lapply(ext3, function(x) bootstats(x,
                                               type = type,
                                               conf = conf)) %>%
      data.table::rbindlist()
  } else {
    Sijk <- NULL
  }
  return(rbind(Si.STi, Sij, Sijk))
}

# CREATE TWO VECTORS WITH SENSITIVITY INDICES AND MODEL INPUTS ----------------

#' Creates two vectors with the sensitivity indices and the model
#' inputs.
#'
#' @param params Vector with the name of the model inputs.
#' @param second Boolean. If second == TRUE, it computes the confidence
#' intervals for second-order indices. Default is second == FALSE.
#' @param third Boolean. If third == TRUE, it computes the confidence
#' intervals for second-order indices. Default is third == FALSE.
#'
#' @return A data table.
create_vectors <- function(params, second = FALSE, third = FALSE) {
  # Define for first and total only
  parameters <- rep(params, each = 2)
  sensitivity <- rep(c("Si", "STi"), times = length(params))
  # Define for second only
  parameters.second <- utils::combn(params, 2, simplify = FALSE) %>%
    lapply(., function(x) paste0(x, collapse = ".")) %>%
    unlist()
  sensitivity.second <- rep("Sij", times = (length(params) * (length(params) - 1) / 2))
  # Define for third only
  if(length(params) > 2) {
    parameters.third <- utils::combn(params, 3, simplify = FALSE) %>%
      lapply(., function(x) paste0(x, collapse = ".")) %>%
      unlist()
    sensitivity.third <- rep("Sijk", times = factorial(length(params)) /
                               (factorial(3) * factorial((length(params) - 3))))
  }
  out <- data.table::as.data.table(cbind(parameters, sensitivity))
  # Set conditions
  if(second == TRUE & third == FALSE) {
    parameters <- c(parameters, parameters.second)
    sensitivity <- c(sensitivity, sensitivity.second)
    out <- data.table::as.data.table(cbind(parameters, sensitivity))
  }
  if(second == FALSE & third == TRUE) {
    parameters <- c(parameters, parameters.third)
    sensitivity <- c(sensitivity, sensitivity.third)
    out <- data.table::as.data.table(cbind(parameters, sensitivity))
  }
  if(second == TRUE & third == TRUE) {
    parameters <- c(parameters, parameters.second, parameters.third)
    sensitivity <- c(sensitivity, sensitivity.second, sensitivity.third)
    out <- data.table::as.data.table(cbind(parameters, sensitivity))
  }
  return(out)
}

#' Computes bootstrap confidence intervals for Sobol' indices
#'
#' @param b A boot object with the computed Sobol' indices.
#' @param params A vector with the name of the model inputs.
#' @param type A vector of character strings representing
#' the type of intervals required. The value should be any subset
#' of the values c('norm','basic', 'perc', 'bca'). For more information,
#' check the function boot::boot.ci.
#' @param conf A scalar or vector containing the confidence
#' level(s) of the required interval(s).
#' @param second Boolean. If second == TRUE, it computes the confidence
#' intervals for second-order indices. Default is second == FALSE.
#' @param third Boolean. If third == TRUE, it computes the confidence
#' intervals for second-order indices. Default is third == FALSE.
#'
#' @return A data table.
#'
#' @export
#'
#' @examples
#' n <- 5000; k <- 8; R <- 10
#' A <- sobol_matrices(n = n, k = k, second = TRUE, third = TRUE)
#' Y <- sobol.Fun(A)
#' sens <- sobol_indices(Y = Y, params = colnames(data.frame(A)),
#' R = R, n = n, parallel = "no", ncpus = 1,
#' second = TRUE, third = TRUE)
#' sobol_ci(sens, params = colnames(data.frame(A)),
#' type = "norm", conf = 0.95, second = TRUE, third = TRUE)
sobol_ci <- function(b, params, type, conf, second = FALSE, third = FALSE) {
  out <- sobol_ci_temp(b = b, params = params, type = type,
                        conf = conf, second = second, third = third) %>%
    cbind(create_vectors(params = params, second = second, third = third))
  return(out)
}
