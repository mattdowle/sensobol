

# PLOT SOBOL' FIRST AND TOTAL-ORDER INDICES -----------------------------------

#' Plot Sobol' first and total-order indices
#'
#' @param x A data.table
#' @param type An integer. If type == 1, it plots first and total effects.
#' If type == 2, it plots second-order effects. If type == 3, it plots
#' third-order effects. Default is type == 1.
#'
#' @return A ggplot object
#' @import ggplot2
#' @export
#'
#' @examples
#' n <- 5000; k <- 8; R <- 10
#' A <- sobol_matrices(n = n, k = k, second = TRUE, third = TRUE)
#' Y <- sobol.Fun(A)
#' sens <- sobol_indices(Y = Y, params = colnames(data.frame(A)),
#' R = R, n = n, parallel = "no", ncpus = 1,
#' second = TRUE, third = TRUE)
#' sens.ci <- sobol_ci(sens, params = colnames(data.frame(A)),
#' type = "norm", conf = 0.95, second = TRUE, third = TRUE) %>%
#' cbind(create.vectors(paste("X", 1:8, sep = ""), second = TRUE, third = TRUE))
#' plot_sobol(sens.ci, type = 2)

plot_sobol <- function(x, type = 1) {
  if(type == 1) {
    gg <- x %>%
      .[sensitivity == "Si" | sensitivity == "STi"] %>%
      ggplot2::ggplot(., aes(parameters, original,
                             fill = sensitivity)) +
      geom_bar(stat = "identity",
               position = position_dodge(0.6),
               color = "black") +
      geom_errorbar(aes(ymin = low.ci,
                        ymax = high.ci),
                    position = position_dodge(0.6)) +
      scale_fill_discrete(name = "Sobol' indices",
                          labels = c(expression(S[italic(i)]),
                                     expression(S[italic(T[i])]))) +
      labs(x = "",
           y = "Variance") +
      theme_bw() +
      theme(legend.position = "top",
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            legend.background = element_rect(fill = "transparent",
                                             color = NA),
            legend.key = element_rect(fill = "transparent",
                                      color = NA))
  } else if(!type == 1) {
    if(type == 2) {
      plot.type <- "Sij"
    } else if(type == 3) {
      plot.type <- "Sijk"
    } else {
      stop("Type should be either 1, 2 or 3")
    }
    gg <- x %>%
      .[sensitivity == plot.type] %>%
      ggplot2::ggplot(., aes(reorder(parameters, original),
                    original)) +
      geom_point() +
      geom_errorbar(aes(ymin = low.ci,
                        ymax = high.ci)) +
      theme_bw() +
      labs(x = "",
           y = "Variance") +
      geom_hline(yintercept = 0,
                 lty = 2,
                 color = "red") +
      theme(panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            legend.background = element_rect(fill = "transparent",
                                             color = NA),
            legend.key = element_rect(fill = "transparent",
                                      color = NA),
            axis.text.x = element_text(angle = 45,
                                       hjust = 1))
  }
  return(gg)
}
