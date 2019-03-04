context("test-sobol_matrices")

N <- 100; k <- 3

test_that("Output is a matrix", {
  A <- sobol_matrices(N, k)
  expect_is(A, "matrix")
})

test_that("Sample size for first and total-order indices", {
  A <- sobol_matrices(N, k)
  expect_equal(nrow(A), N * (k + 2))
})

test_that("Sample size for first, total and second-order indices", {
  A <- sobol_matrices(N, k, second = TRUE)
  expect_equal(nrow(A), (N * (k + 2)) + (N * factorial(k) /
                                           factorial(2) * factorial(k - 2)))
})

test_that("Sample size for first, total, second and third-order indices", {
  A <- sobol_matrices(N, k, second = TRUE, third = TRUE)
  expect_equal(nrow(A), (N * (k + 2)) + (N * factorial(k) /
                                           factorial(2) * factorial(k - 2)) +
                 (N * factorial(k) /
                    factorial(3) * factorial(k - 3)))
})
