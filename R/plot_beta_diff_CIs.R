#' Plot credible intervals for adjacent coefficient differences
#'
#' @param Beta Draws-by-covariates-by-classes array including the reference
#'   class, or a fit from [bcmlr()].
#' @param coef_names Optional distinct covariate labels in input-column order.
#' @param ci_level Equal-tailed credible-interval level in `(0, 1)`.
#' @param ncol Number of facet columns; defaults to the number of changepoints.
#' @param title Optional plot title.
#' @return A ggplot object. Its data contain posterior means, interval endpoints,
#'   covariate labels and the ordered adjacent contrast index.
#' @details Differences are computed draw by draw as `beta[, j+1] - beta[, j]`.
#'   Coefficients refer to the standardized covariates used for fitting.
#' @export
plot_beta_diff_CIs <- function(Beta, coef_names = NULL, ci_level = 0.95,
                               ncol = NULL, title = NULL) {
  defaults <- if (is.list(Beta)) colnames(Beta$data) else dimnames(Beta)[[2L]]
  B <- .beta_draws(Beta)
  p <- dim(B)[2L]
  L <- dim(B)[3L] - 1L
  coef_names <- .labels(coef_names, p, defaults)
  .scalar(ci_level, "ci_level", 0, 1, open_lower = TRUE, open_upper = TRUE)
  if (is.null(ncol)) ncol <- L
  .scalar(ncol, "ncol", 1, integer = TRUE)
  labs <- paste0("beta[", seq_len(L) + 1L, "] - beta[", seq_len(L), "]")
  d <- do.call(rbind, lapply(seq_len(L), function(j) {
    delta <- matrix(B[, , j + 1L] - B[, , j], nrow = dim(B)[1L], ncol = p)
    ci <- vapply(seq_len(p), function(k) stats::quantile(delta[, k],
                 c((1 - ci_level) / 2, (1 + ci_level) / 2), names = FALSE), numeric(2))
    data.frame(coef = coef_names, mean = colMeans(delta), lower = ci[1, ], upper = ci[2, ],
                contrast = j, label = labs[j])
  }))
  d$coef <- factor(d$coef, levels = coef_names)
  d$label <- factor(d$label, levels = labs)
  ggplot2::ggplot(d, ggplot2::aes(x = .data$mean, y = .data$coef)) +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed") +
    ggplot2::geom_segment(ggplot2::aes(x = .data$lower, xend = .data$upper, yend = .data$coef)) +
    ggplot2::geom_point(size = 1.8) +
    ggplot2::facet_wrap(~ label, ncol = ncol, labeller = ggplot2::label_parsed) +
    ggplot2::labs(x = NULL, y = NULL, title = title) + ggplot2::theme_bw()
}
