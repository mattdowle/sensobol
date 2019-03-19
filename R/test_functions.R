
# SOBOL' G function -----------------------------------------------------------

#' Sobol' G function
#'
#' @param X A data frame or matrix.
#'
#' @return A numeric vector with the model output.
#' @export
#'
#' @examples
#' A <- sobol_matrices(n = 100, k = 8)
#' Y <- sobol_Fun(A)
sobol_Fun <- function(X) {
  a <- c(0, 1, 4.5, 9, 99, 99, 99, 99)
  y <- 1
  for (j in 1:8) {
    y <- y * (abs(4 * X[, j] - 2) + a[j])/(1 + a[j])
  }
  return(y)
}

# Ishigami function -----------------------------------------------------------

ishigami <- function(X1, X2, X3) {
  A <- 2
  B <- 1
  sin(X1) + A * sin(X2) ^ 2 + B * X3 ^ 4 * sin(X1)
}

#' Ishigami function
#'
#' @param X A data frame, data table or matrix with the three model inputs
#' required to run the Ishigami function.
#'
#' @return A numeric vector with the model output.
#' @export
#'
#' @examples
#' A <- sobol_matrices(n = 100, k = 3)
#' Y <- ishigami_Mapply(A)
ishigami_Mapply <- function(X) {
  return(mapply(ishigami,
                X[, 1],
                X[, 2],
                X[, 3]))
}

# Bradley et al. function -----------------------------------------------------


#' Bratley et al. (1992) function
#'
#' It implements the \insertCite{Bratley1992;textual}{sensobol} function.
#'
#' @param X A data frame or numeric matrix.
#'
#' @return A numeric vector with the model output.
#' @export
#' @references
#' \insertAllCited{}
#'
#' @examples
#' A <- sobol_matrices(n = 100, k = 4)
#' Y <- bratley_Fun(A)
bratley_Fun <- function(X) {
  # Preallocate
  xxmat <- xxmatlow <- tmp <- vector(mode = "list",
                                     length = nrow(X))
  Y <- vector(mode = "numeric",
              length = nrow(X))
  for(i in 1:nrow(X)) {
    xxmat[[i]] <- matrix(rep(X[i, ], times = ncol(X)),
                         nrow = ncol(X),
                         ncol = ncol(X),
                         byrow = TRUE)
    xxmatlow[[i]] <- xxmat[[i]]
    xxmatlow[[i]][upper.tri(xxmat[[i]])] <- 1
    tmp[[i]] <- matrixStats::rowProds(xxmatlow[[i]])
    Y[[i]] <- sum(tmp[[i]] * (-1) ^ (1:ncol(X)))
  }
  return(Y)
}
