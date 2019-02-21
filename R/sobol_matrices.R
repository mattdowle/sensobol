
# INTERNAL FUNCTION TO CREATE THE SCRAMBLED MATRIX FOR THE
# COMPUTATION OF FIRST-ORDER SOBOL' INDICES -----------------------------------

#' Creation of Sobol' matrices for the computation of
#' first and total order effects. It is an internal function of
#' "sobol_matrices".
#'
#' @param A The first k Sobol' matrix.
#' @param B The second k Sobol' matrix.
#'
#' @return A matrix.
#' @export

scrambled_sobol <- function(A, B) {
  X <- rbind(A, B)
  for(i in 1:ncol(A)) {
    AB <- A
    AB[, i] <- B[, i]
    X <- rbind(X, AB)
  }
  AB <- X
  return(AB)
}

# INTERNAL FUNCTION TO CREATE THE SCRAMBLED MATRIX FOR THE
# COMPUTATION OF SECOND-ORDER SOBOL' INDICES ----------------------------------

#' Creation of Sobol' matrix for the computation of third-order effects. It is
#' an internal function of "sobol_matrices".
#'
#' @param A First k Sobol' matrix.
#' @param B Second k Sobol' matrix.
#'
#' @return A matrix
#' @export

scrambled_second <- function(A, B) {
  X <- rbind(A, B)
  parms <- utils::combn(1:ncol(A), 2, simplify = FALSE)
  for(i in parms) {
    AB <- A
    AB[, c(i[[1]], i[[2]])] <- B[, c(i[[1]], i[[2]])]
    X <- rbind(X, AB)
  }
  AB <- X %>%
    .[(2 * nrow(A) + 1):nrow(.), ]
  return(AB)
}

# INTERNAL FUNCTION TO CREATE THE SCRAMBLED MATRIX FOR THE
# COMPUTATION OF THIRD-ORDER SOBOL' INDICES ----------------------------------

#' Creation of Sobol' matrix for the computation of third-order effects. It is
#' an internal function of "sobol_matrices".
#'
#' @param A First k Sobol' matrix.
#' @param B Second k Sobol' matrix.
#'
#' @return A matrix
#' @export
#'
scrambled_third <- function(A, B) {
  X <- rbind(A, B)
  parms <- utils::combn(1:ncol(A), 3, simplify = FALSE)
  for(i in parms) {
    AB <- A
    AB[, c(i[[1]], i[[2]], i[[3]])] <- B[, c(i[[1]], i[[2]], i[[3]])]
    X <- rbind(X, AB)
  }
  AB <- X %>%
    # Select only the scrambled matrix
    .[(2 * nrow(A) + 1):nrow(.), ]
  return(AB)
}

# FUNCTION TO CREATE THE SOBOL' MATRIX FOR THE COMPUTATION OF
# FIRST AND TOTAL ORDER EFFECTS - AS WELL AS SECOND AND THIRD
# ORDER EFFECTS, IF DESIRED ---------------------------------------------------

#' Creation of Sobol' matrices to compute first, second, third and total-
#' order effects.
#'
#' @param n Integer, sample size of the Sobol' matrix.
#' @param k Integer, number of model inputs.
#' @param second Boolean. If second == TRUE, it creates the scrambled
#' matrix required to compute second-order effects. Default is second == FALSE.
#' @param third Boolean. If third == TRUE, it creates the scrambled
#' matrix required to compute third-order effects. Default is second == FALSE.
#'
#' @return A matrix
#' @export
#'
#' @examples
#' sobol_matrices(n = 1000, k = 8, second = TRUE, third = TRUE)
sobol_matrices <- function(n, k, second = FALSE, third = FALSE) {
  # Create the Sobol quasi-random number sequence
  df <- randtoolbox::sobol(n = n,
                           dim = k * 2)
  # Create A matrix
  A <- df[, 1:k]
  # Create B matrix
  B <- df[, (k + 1) : (k * 2)]
  # Create AB matrix
  AB <- scrambled_sobol(A = A, B = B)
  if(second == TRUE & third == FALSE) {
    # Compute AB matrix for second order indices
    AB.2 <- scrambled_second(A = A, B = B)
    AB <- rbind(AB, AB.2)
  }
  if(second == FALSE & third == TRUE) {
    AB.3 <- scrambled_third(A = A, B = B)
    AB <- rbind(AB, AB.3)
  }
  if(second == TRUE & third == TRUE) {
    AB.2 <- scrambled_second(A = A, B = B)
    AB.3 <- scrambled_third(A = A, B = B)
    AB <- rbind(AB, AB.2, AB.3)
  }
  return(AB)
}
