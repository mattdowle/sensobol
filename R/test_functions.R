
# SOBOL' G function -----------------------------------------------------------

#' Sobol' G function
#'
#' @param X A data frame or matrix.
#'
#' @return A numeric vector with the model output.
#' @export
#'
#' @examples
#' A <- sobol_matrices(n = 1000, k = 8)
#' Y <- sobol.Fun(A)
sobol.Fun <- function(X) {
  a <- c(0, 1, 4.5, 9, 99, 99, 99, 99)
  y <- 1
  for (j in 1:8) {
    y <- y * (abs(4 * X[, j] - 2) + a[j])/(1 + a[j])
  }
  return(y)
}
