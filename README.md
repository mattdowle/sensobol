# sensobol

The goal of `sensobol` is to provide a set of functions to swiftly compute and visualize Sobol' sensitivity indices. The functions allow to: 
- Create the sample matrices for the model evaluation.
- Swiftly compute and bootstrap up to third-order effects.
- Assess the approximation error of Sobol' indices (Sobol' indices of a dummy parameter)
- Plot the model uncertainty and the Sobol' indices.

## Installation

You can install the released version of sensobol as follows:

``` r
install.packages("devtools") # if you have not installed devtools package already
devtools::install_github("arnaldpuy/sensobol", build_vignettes = TRUE)
```

## Example

``` r
## Create sample matrix
A <- sobol_matrices(n = 1000, k = 3, second = TRUE)

## Compute the model output:
Y <- ishigami_Mapply(A)

## Compute the Sobol' indices:
sens <- sobol_indices(Y = Y, params = colnames(data.frame(A)),
R = 100, n = 1000, second = TRUE)

## Plot Sobol' indices:
plot_sobol(sens.ci, type = 1)
```
