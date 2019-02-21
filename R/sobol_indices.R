

# FUNCTION TO COMPUTE SOBOL' FIRST AND TOTAL-ORDER EFFECTS USING
# THE JANSEN 1999 ESTIMATOR FOR FIRST AND TOTAL INDICES -----------------------

#' Computation of first and total-effect Sobol' indices using the
#' Jansen 1999 estimators.
#'
#' @param Y_A Numeric vector, model output of the A matrix.
#' @param Y_B Numeric vector, model output of the B matrix.
#' @param Y_AB Numeric vector, model output of the AB matrix,
#' where the AB matrix is composed of all columns of the
#' A matrix except the i-th, which comes from the B matrix.
#'
#' @return A numeric vector with the first and total Sobol' indices.
sobol_computeJ <- function(Y_A, Y_B, Y_AB) {
  # Compute sample mean of output
  f0 <- (1 / (2 * length(Y_A))) * sum(Y_A + Y_B)
  # Compute unconditional variance of output
  VY <- 1 / (2 * length(Y_A) - 1) * sum((Y_A - f0) ^ 2 + (Y_B - f0) ^ 2)
  # Compute first order indices (Jansen 1999)
  Si <- (VY - ((1 / (2 * length(Y_A))) * sum((Y_B - Y_AB) ^ 2))) / VY
  # Compute total order indices (Jansen 1999)
  STi <- ((1 / (2 * length(Y_A))) * sum((Y_A - Y_AB) ^ 2)) / VY
  return(c(Si, STi))
}

#' Wrapper to allow bootstrap of the sobol_computeJ function
#'
#' @param d Data table with the model output of
#' the A, B and AB matrix as columns
#' @param i Indices to allow bootstrap.
#'
#' @return The first and total Sobol' indices.
sobol_MapplyJ <- function(d, i) {
  return(mapply(sobol_computeJ,
                d[i, "Y_A"],
                d[i, "Y_B"],
                d[i, "Y_AB"]))
}

# FUNCTION TO COMPUTE SOBOL' SECOND-ORDER EFFECTS USING
# THE JANSEN 1999 ESTIMATOR ---------------------------------------------------

#' Computation of second-order Sobol' indices using the
#' Jansen 1999 estimators.
#'
#' @param Y_A Numeric vector, model output of the A matrix.
#' @param Y_B Numeric vector, model output of the B matrix.
#' @param Y_ABi Numeric vector, model output of the AB matrix,
#' where all the columns come from the A matrix except the i-th,
#' which comes from the B matrix.
#' @param Y_ABj Numeric vector, model output of the AB matrix,
#' where all the columns come from the A matrix except the j-th,
#' which comes from the B matrix.
#' @param Y_ABij Numeric vector, model output of the AB matrix,
#' where all the columns come from the A matrix except the i-th
#' and the j-th, which come from the B matrix.
#'
#' @return The second-order indices
sobol_compute_secondJ <- function(Y_A, Y_B, Y_ABi,
                                  Y_ABj, Y_ABij) {
  # Compute sample mean of output
  f0 <- (1 / (2 * length(Y_A))) * sum(Y_A + Y_B)
  # Compute unconditional variance of output
  VY <- 1 / (2 * length(Y_A) - 1) * sum((Y_A - f0) ^ 2 + (Y_B - f0) ^ 2)
  # Use Jansen 1999 estimates
  Vi <- VY - ((1 / (2 * length(Y_A))) * sum((Y_B - Y_ABi) ^ 2))
  Vj <- VY - ((1 / (2 * length(Y_A))) * sum((Y_B - Y_ABj) ^ 2))
  Vij <- (VY - ((1 / (2 * length(Y_A))) * sum((Y_B - Y_ABij) ^ 2))) - Vi - Vj
  Sij <- Vij / VY
  return(Sij)
}

# Wrapper to allow bootstrap
#' Wrapper to allow bootstrap of the sobol_compute_secondJ function
#'
#' @param d Data table with the model output of
#' the A, B and AB matrix as columns
#' @param i Indices to allow bootstrap.
#'
#' @return The second-order indices.
sobol_second_MapplyJ <- function(d, i) {
  return(mapply(sobol_compute_secondJ,
                d[i, "Y_A"],
                d[i, "Y_B"],
                d[i, "Y_ABi"],
                d[i, "Y_ABj"],
                d[i, "Y_ABij"]))
}

# FUNCTION TO COMPUTE SOBOL' THIRD-ORDER EFFECTS USING
# THE JANSEN 1999 ESTIMATOR ---------------------------------------------------

#' Computation of third-order Sobol' indices using the
#' Jansen 1999 estimators.
#'
#' @param Y_A Numeric vector, model output of the A matrix.
#' @param Y_B Numeric vector, model output of the B matrix.
#' @param Y_ABi Numeric vector, model output of the AB matrix,
#' where all the columns come from the A matrix except the i-th,
#' which comes from the B matrix.
#' @param Y_ABj Numeric vector, model output of the AB matrix,
#' where all the columns come from the A matrix except the j-th,
#' which comes from the B matrix.
#' @param Y_ABk Numeric vector, model output of the AB matrix,
#' where all the columns come from the A matrix except the k-th,
#' which comes from the B matrix.
#' @param Y_ABij Numeric vector, model output of the AB matrix,
#' where all the columns come from the A matrix except the i-th
#' and the j-th, which come from the B matrix.
#' @param Y_ABik Numeric vector, model output of the AB matrix,
#' where all the columns come from the A matrix except the i-th
#' and the k-th, which come from the B matrix.
#' @param Y_ABjk Numeric vector, model output of the AB matrix,
#' where all the columns come from the A matrix except the j-th
#' and the k-th, which come from the B matrix.
#' @param Y_ABijk Numeric vector, model output of the AB matrix,
#' where all the columns come from the A matrix except the i-th,
#' the j-th and the k-th, which come from the B matrix.
#'
#' @importFrom stats "na.omit"
#' @return The third-order indices.
sobol_compute_thirdJ <- function(Y_A, Y_B, Y_ABi,
                                 Y_ABj, Y_ABk, Y_ABij,
                                 Y_ABik, Y_ABjk, Y_ABijk) {
  # Compute sample mean of output
  f0 <- (1 / (2 * length(Y_A))) * sum(Y_A + Y_B)
  # Compute unconditional variance of output
  VY <- 1 / (2 * length(Y_A) - 1) * sum((Y_A - f0) ^ 2 + (Y_B - f0) ^ 2)
  # Use Jansen 1999 estimates
  Vi <- VY - ((1 / (2 * length(Y_A))) * sum((Y_B - Y_ABi) ^ 2))
  Vj <- VY - ((1 / (2 * length(Y_A))) * sum((Y_B - Y_ABj) ^ 2))
  Vk <- VY - ((1 / (2 * length(Y_A))) * sum((Y_B - Y_ABk) ^ 2))
  Vij <- (VY - ((1 / (2 * length(Y_A))) * sum((Y_B - Y_ABij) ^ 2))) - Vi - Vj
  Vik <- (VY - ((1 / (2 * length(Y_A))) * sum((Y_B - Y_ABik) ^ 2))) - Vi - Vk
  Vjk <- (VY - ((1 / (2 * length(Y_A))) * sum((Y_B - Y_ABjk) ^ 2))) - Vj - Vk
  Vijk <- (VY - ((1 / (2 * length(Y_A))) * sum((Y_B - Y_ABijk) ^ 2))) - Vij - Vik - Vjk - Vi - Vj - Vk
  Sijk <- Vijk / VY
  return(Sijk)
}

# Wrap-up the Sobol' function for third order
# indices and allow bootstrap
#' Wrapper to allow bootstrap of the sobol_compute_thirdJ function
#'
#' @param d Data table with the model output of
#' the A, B and AB matrix as columns
#' @param i Indices to allow bootstrap.
#'
#' @return The third-order indices.
sobol_third_MapplyJ <- function(d, i) {
  return(mapply(sobol_compute_thirdJ,
                d[i, "Y_A"],
                d[i, "Y_B"],
                d[i, "Y_ABi"],
                d[i, "Y_ABj"],
                d[i, "Y_ABk"],
                d[i, "Y_ABij"],
                d[i, "Y_ABik"],
                d[i, "Y_ABjk"],
                d[i, "Y_ABijk"]))
}

# FUNCTION TO COMPUTE SOBOL' FIRST AND TOTAL-ORDER EFFECTS USING
# THE JANSEN 1999 ESTIMATOR FOR FIRST AND TOTAL INDICES -----------------------

#' Computation of first and total-effect Sobol' indices using the
#' Saltelli et al. 2010 estimators.
#'
#' @param Y_A Numeric vector, model output of the A matrix.
#' @param Y_B Numeric vector, model output of the B matrix.
#' @param Y_AB Numeric vector, model output of the AB matrix,
#' where the AB matrix is composed of all columns of the
#' A matrix except the i-th, which comes from the B matrix.
#'
#' @return A numeric vector with the first and total Sobol' indices.

sobol_computeS <- function(Y_A, Y_B, Y_AB) {
  # Compute sample mean of output
  f0 <- (1 / (2 * length(Y_A))) * sum(Y_A + Y_B)
  # Compute unconditional variance of output
  VY <- 1 / (2 * length(Y_A) - 1) * sum((Y_A - f0) ^ 2 + (Y_B - f0) ^ 2)
  # Compute first order indices (Saltelli et al. 2010)
  Si <- (1 / length(Y_A)) * sum(Y_B * (Y_AB - Y_A)) / VY
  # Compute total order indices (Jansen 1999)
  STi <- ((1 / (2 * length(Y_A))) * sum((Y_A - Y_AB) ^ 2)) / VY
  return(c(Si, STi))
}



#' Wrapper to allow bootstrap of the sobol_computeS function
#'
#' @param d Data table with the model output of
#' the A, B and AB matrix as columns
#' @param i Indices to allow bootstrap.
#'
#' @return The first and total Sobol' indices.
sobol_MapplyS <- function(d, i)
  return(mapply(sobol_computeS,
                d[i, "Y_A"],
                d[i, "Y_B"],
                d[i, "Y_AB"]))

# FUNCTION TO COMPUTE SOBOL' SECOND-ORDER EFFECTS USING
# THE SALTELLI ET AL. 2010 ESTIMATOR ------------------------------------------

#' Computation of second-order Sobol' indices using the
#' Saltelli et al. 2010 estimators.
#'
#' @param Y_A Numeric vector, model output of the A matrix.
#' @param Y_B Numeric vector, model output of the B matrix.
#' @param Y_ABi Numeric vector, model output of the AB matrix,
#' where all the columns come from the A matrix except the i-th,
#' which comes from the B matrix.
#' @param Y_ABj Numeric vector, model output of the AB matrix,
#' where all the columns come from the A matrix except the j-th,
#' which comes from the B matrix.
#' @param Y_ABij Numeric vector, model output of the AB matrix,
#' where all the columns come from the A matrix except the i-th
#' and the j-th, which come from the B matrix.
#'
#' @return The second-order indices
sobol_compute_secondS <- function(Y_A, Y_B, Y_ABi,
                                  Y_ABj, Y_ABij) {
  # Compute sample mean of output
  f0 <- (1 / (2 * length(Y_A))) * sum(Y_A + Y_B)
  # Compute unconditional variance of output
  VY <- 1 / (2 * length(Y_A) - 1) * sum((Y_A - f0) ^ 2 + (Y_B - f0) ^ 2)
  # Use Saltelli et al. 2010 estimate
  Vi <- (1 / length(Y_A)) * sum(Y_B * (Y_ABi - Y_A))
  Vj <- (1 / length(Y_A)) * sum(Y_B * (Y_ABj - Y_A))
  Vij <- ((1 / length(Y_A)) * sum(Y_B * (Y_ABij - Y_A))) - Vi - Vj
  Sij <- Vij / VY
  return(Sij)
}


#' Wrapper to allow bootstrap of the sobol_compute_secondS function
#'
#' @param d Data table with the model output of
#' the A, B and AB matrix as columns
#' @param i Indices to allow bootstrap.
#'
#' @return The second-order indices.
sobol_second_MapplyS <- function(d, i)
  return(mapply(sobol_compute_secondS,
                d[i, "Y_A"],
                d[i, "Y_B"],
                d[i, "Y_ABi"],
                d[i, "Y_ABj"],
                d[i, "Y_ABij"]))

# FUNCTION TO COMPUTE SOBOL' THIRD-ORDER EFFECTS USING
# THE SALTELLI ET AL. 2010 ESTIMATOR ------------------------------------------

#' Computation of third-order Sobol' indices using the
#' Saltelli et al. 2010 estimators.
#'
#' @param Y_A Numeric vector, model output of the A matrix.
#' @param Y_B Numeric vector, model output of the B matrix.
#' @param Y_ABi Numeric vector, model output of the AB matrix,
#' where all the columns come from the A matrix except the i-th,
#' which comes from the B matrix.
#' @param Y_ABj Numeric vector, model output of the AB matrix,
#' where all the columns come from the A matrix except the j-th,
#' which comes from the B matrix.
#' @param Y_ABk Numeric vector, model output of the AB matrix,
#' where all the columns come from the A matrix except the k-th,
#' which comes from the B matrix.
#' @param Y_ABij Numeric vector, model output of the AB matrix,
#' where all the columns come from the A matrix except the i-th
#' and the j-th, which come from the B matrix.
#' @param Y_ABik Numeric vector, model output of the AB matrix,
#' where all the columns come from the A matrix except the i-th
#' and the k-th, which come from the B matrix.
#' @param Y_ABjk Numeric vector, model output of the AB matrix,
#' where all the columns come from the A matrix except the j-th
#' and the k-th, which come from the B matrix.
#' @param Y_ABijk Numeric vector, model output of the AB matrix,
#' where all the columns come from the A matrix except the i-th,
#' the j-th and the k-th, which come from the B matrix.
#'
#' @importFrom stats "na.omit"
#' @return The third-order indices.
sobol_compute_thirdS <- function(Y_A, Y_B, Y_ABi,
                                 Y_ABj, Y_ABk, Y_ABij,
                                 Y_ABik, Y_ABjk, Y_ABijk) {
  # Compute sample mean of output
  f0 <- (1 / (2 * length(Y_A))) * sum(Y_A + Y_B)
  # Compute unconditional variance of output
  VY <- 1 / (2 * length(Y_A) - 1) * sum((Y_A - f0) ^ 2 + (Y_B - f0) ^ 2)
  # Use Saltelli et al. 2010 estimate
  Vi <- (1 / length(Y_A)) * sum(Y_B * (Y_ABi - Y_A))
  Vj <- (1 / length(Y_A)) * sum(Y_B * (Y_ABj - Y_A))
  Vk <- (1 / length(Y_A)) * sum(Y_B * (Y_ABk - Y_A))
  Vij <- ((1 / length(Y_A)) * sum(Y_B * (Y_ABij - Y_A))) -  Vi - Vj
  Vik <- ((1 / length(Y_A)) * sum(Y_B * (Y_ABik - Y_A))) - Vi - Vk
  Vjk <- ((1 / length(Y_A)) * sum(Y_B * (Y_ABjk - Y_A))) - Vj - Vk
  Vijk <- ((1 / length(Y_A)) * sum(Y_B * (Y_ABijk - Y_A))) - Vij - Vik - Vjk - Vi - Vj - Vk
  Sijk <- Vijk / VY
  return(Sijk)
}

#' Wrapper to allow bootstrap of the sobol_compute_thirdS function
#'
#' @param d Data table with the model output of
#' the A, B and AB matrix as columns
#' @param i Indices to allow bootstrap.
#'
#' @return The third-order indices.
# Wrap-up the Sobol' function for second order
# indices and allow bootstrap
sobol_third_MapplyS <- function(d, i)
  return(mapply(sobol_compute_thirdS,
                d[i, "Y_A"],
                d[i, "Y_B"],
                d[i, "Y_ABi"],
                d[i, "Y_ABj"],
                d[i, "Y_ABk"],
                d[i, "Y_ABij"],
                d[i, "Y_ABik"],
                d[i, "Y_ABjk"],
                d[i, "Y_ABijk"]))

# FUNCTION TO COMPUTE FIRST, SECOND, THIRD AND TOTAL
# SOBOL' INDICES --------------------------------------------------------------

#' Computation of first, second, third and total-effects Sobol' indices.
#'
#' It computes and bootstraps up to third-order Sobol' indices
#' using either the Saltelli et al. 2010 or the Jansen 1999
#' estimator. It relies on parallel computing.
#'
#' @param Y Numeric vector, model output.
#' @param params Vector with the name of the model inputs.
#' @param type Estimator to use. type == "saltelli" to use
#' the Saltelli et al. 2010 estimator, or type == "jansen"
#' to use the Jansen 1999 estimator. Default is type == "jansen.
#' @param R Integer, number of bootstrap replicas.
#' @param n Integer, sample size of the sample matrix.
#' @param parallel The type of parallel operation to be used (if any).
#' If missing, the default is taken from the option "boot.parallel"
#' (and if that is not set, "no"). For more information, check the
#' "parallel" option in the boot::boot function.
#' @param ncpus Integer: number of processes to be used in parallel operation:
#' typically one would chose this to the number of available CPUs.
#' For more information, check the "parallel" option in
#' the boot::boot function.
#' @param second Boolean. if second == "TRUE", it computes
#' second-order Sobol' indices.
#' @param third Boolean. if third == "TRUE", it computes
#' third-order Sobol' indices.
#' @importFrom data.table ".SD"
#'
#' @return A data.table object.
#' @export
#' @examples
#' n <- 5000; k <- 8; R <- 10
#' A <- sobol_matrices(n = n, k = k, second = TRUE, third = TRUE)
#' Y <- sobol.Fun(A)
#' sens <- sobol_indices(Y = Y, params = colnames(data.frame(A)),
#' R = R, n = n, parallel = "no", ncpus = 1,
#' second = TRUE, third = TRUE)
sobol_indices <- function(Y, params, type = "jansen",
                          R, n, parallel = "no", ncpus = 1,
                          second = FALSE, third = FALSE) {
  # Calculate the number of parameters
  k <- length(params)
  # Calculate the length of the A and B matrices
  p <- length(1:n)
  # Extract the model output of the A matrix
  Y_A <- Y[1:p]
  # Extract the model output of the B matrix
  Y_B <- Y[(p + 1) : (2 * p)]
  # Extract the model output of the AB matrix
  Y_AB <- Y[(2*p+1):(n * (k + 2))]
  # Create vector with parameters
  parameters <- rep(params, each = length(Y_A))
  # merge vector with data table
  vec <- cbind(Y_A, Y_B, Y_AB)
  out <- data.table::data.table(vec, parameters)
  out.1 <- out %>%
    # remove rows with NA
    na.omit(.)
  # Check which estimator to use for First-order indices
  if(type == "jansen") {
    Estimator1 <- sobol_MapplyJ
    Estimator2 <- sobol_second_MapplyJ
    Estimator3 <- sobol_third_MapplyJ
  } else if(type == "saltelli") {
    Estimator1 <- sobol_MapplyS
    Estimator2 <- sobol_second_MapplyS
    Estimator3 <- sobol_third_MapplyS
  }
  # Bootstrap Sobol'indices
  Si.STi <- out.1[, list(list(boot::boot(.SD,
                                         Estimator1,
                                         R = R,
                                         parallel = parallel,
                                         ncpus = ncpus))),
                  by = parameters]
  # CHECK IF THE MODEL OUTPUT INCLUDES
  # MODEL OUTPUT FOR THE CALCULATION OF
  # SECOND ORDER INDICES
  if(second == TRUE) {
    # Extract the model output of the
    # second order matrix
    Y_ABij <- Y[((n * (k + 2)) + 1):
                  ((n * (k + 2)) +
                     (n * factorial(k) / (factorial(2) * factorial(k - 2))))]
    # Create vector with pairs of parameters
    paramet <- utils::combn(params, 2, simplify = FALSE) %>%
      lapply(., function(x) paste0(x, collapse = ".")) %>%
      unlist()
    parameters <- rep(paramet, each = length(Y_A))
    # create vectpr with Y_A, Y_B and Y_ABij
    vec <- cbind(Y_A, Y_B, Y_ABij)
    # Create data table
    out2 <- data.table::data.table(vec, parameters)
    # Set key in the matrix of the first and total indices
    data.table::setkey(out, "parameters")
    # Extract the first parameter
    first <- sub("\\..*", "", paramet)
    # Extract the second parameter
    second <- sub(".*\\.", "", paramet)
    # Extract the AB vector of the first parameter
    Y_ABi <- out[.(first), allow.cartesian = TRUE] %>%
      .[, Y_AB]
    # Extract the AB vector of the second parameter
    Y_ABj <- out[.(second), allow.cartesian = TRUE] %>%
      .[, Y_AB]
    # Merge with the AB matrix of the second order
    out2[, Y_ABi:= cbind(Y_ABi)] %>%
      .[, Y_ABj:= cbind(Y_ABj)]
    # Remove rows with NA
    out3 <- out2 %>%
      na.omit(.)
    # Bootstrap second-order indices
    Sij <- out3[, list(list(boot::boot(.SD,
                                       Estimator2,
                                       R = R,
                                       parallel = parallel,
                                       ncpus = ncpus))),
                by = parameters]
  } else {
    Sij <- NULL
  }
  if(third == TRUE) {
    # Extract the model output of the
    # third order matrix
    Y_ABijk <- Y[(((n * (k + 2)) + (n * factorial(k) /
                                      (factorial(2) * factorial(k - 2))))+1):
                   length(Y)]
    # Create vector with pairs of parameters
    paramet3 <- utils::combn(params, 3, simplify = FALSE) %>%
      lapply(., function(x) paste0(x, collapse = ".")) %>%
      unlist()
    parameters <- rep(paramet3, each = length(Y_A))
    # create vectpr with Y_A, Y_B and Y_ABijk
    vec <- cbind(Y_A, Y_B, Y_ABijk)
    # Create data table
    out4 <- data.table::data.table(vec, parameters)
    # Set key in the matrix of the first and total indices
    data.table::setkey(out, "parameters")
    # Set key in the matrix of the second-order indices
    data.table::setkey(out2, "parameters")
    # Extract the first parameter
    first <- sub("\\..*", "", paramet3)
    # Extract the second parameter
    second <- unlist(lapply(paramet3, function(x) unlist(strsplit(x, "[.]"))[[2]]))
    # Extract the third parameter
    third <- sub(".*\\.", "", paramet3)
    # Extract the ij parameter
    ij <- sub(".[^.]*$", "", paramet3)
    # Extract the jk parameter
    jk <- gsub("^.*?\\.","", paramet3)
    # Extract the ik parameter
    ik <- paste(stringr::word(paramet3, 1, sep = stringr::fixed(".")),
                stringr::word(paramet3, 3, sep = stringr::fixed(".")),
                sep = ".")
    # Extract the AB vector of the first parameter
    Y_ABi <- out[.(first), allow.cartesian = TRUE] %>%
      .[, Y_AB]
    # Extract the AB vector of the second parameter
    Y_ABj <- out[.(second), allow.cartesian = TRUE] %>%
      .[, Y_AB]
    # Extract the AB vector of the third parameter
    Y_ABk <- out[.(third), allow.cartesian = TRUE] %>%
      .[, Y_AB]
    # Extract the AB vector of the ij
    Y_ABij <- out2[.(ij), allow.cartesian = TRUE] %>%
      .[, Y_ABij]
    # Extract the AB vector of the jk
    Y_ABjk <- out2[.(jk), allow.cartesian = TRUE] %>%
      .[, Y_ABij]
    # Extract the AB vector of the ik
    Y_ABik <- out2[.(ik), allow.cartesian = TRUE] %>%
      .[, Y_ABij]
    # Merge with the AB matrix of the second order
    out4[, Y_ABi:= cbind(Y_ABi)] %>%
      .[, Y_ABj:= cbind(Y_ABj)] %>%
      .[, Y_ABk:= cbind(Y_ABk)] %>%
      .[, Y_ABij:= cbind(Y_ABij)] %>%
      .[, Y_ABjk:= cbind(Y_ABjk)] %>%
      .[, Y_ABik:= cbind(Y_ABik)]
    # Remove rows with NA
    out5 <- out4 %>%
      na.omit(.)
    # Bootstrap third-order indices
    Sijk <- out5[, list(list(boot::boot(.SD,
                                        Estimator3,
                                        R = R,
                                        parallel = parallel,
                                        ncpus = ncpus))),
                 by = parameters]
  } else {
    Sijk <- NULL
  }
  return(rbind(Si.STi, Sij, Sijk))
}


