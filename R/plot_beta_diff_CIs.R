#' Plot posterior intervals for adjacent coefficient differences
#'
#' @param Beta Numeric array of posterior draws by network statistics by
#'   chronological segments, as returned in `fit$Beta`.
#' @param mfrow Panel layout, c(rows, columns). The paper default is c(1, 4).
#' @param coef_labels Optional labels in the same order as the columns of `S`.
#' @param ci_level Credible-interval level; default 0.95.
#' @param xlim Optional x-axis limits. As in the paper script, these are scale
#'   limits, so out-of-range graphical elements can be removed by ggplot2.
#' @param base_size Base text size for `theme_bw()`.
#' @returns A ggplot object. Points show posterior means of draw-wise
#'   next-minus-previous coefficient differences; bars show equal-tailed intervals.
#' @details The original panel layout, font settings, and axis labels are retained.
#' @examples
#' Beta <- array(0, c(5, 2, 2))
#' Beta[, 1, 2] <- seq(0.1, 0.5, length.out = 5)
#' plot_beta_diff_CIs(Beta, mfrow = c(1, 1), coef_labels = c("Edges", "Triangles"))
#' @export
plot_beta_diff_CIs <- function(Beta, mfrow = c(1, 4), coef_labels = NULL,
                               ci_level = 0.95, xlim = NULL,
                               base_size = 12) {
  # Beta: posterior draws x coefficients x classes
  if (length(dim(Beta)) != 3) stop("Beta must be a 3D array: draws x p x J.")

  p <- dim(Beta)[2]
  J <- dim(Beta)[3]
  if (J < 2) stop("Beta must have at least two classes.")
  if (is.null(coef_labels)) coef_labels <- paste0("coef ", seq_len(p))
  if (length(coef_labels) != p) stop("coef_labels must have length dim(Beta)[2].")

  # Consecutive differences: beta_2 - beta_1, beta_3 - beta_2, ...
  beta_diff <- Beta[, , -1, drop = FALSE] - Beta[, , -J, drop = FALSE]

  probs <- c((1 - ci_level) / 2, 1 - (1 - ci_level) / 2)
  means <- apply(beta_diff, c(2, 3), mean, na.rm = TRUE)
  lower <- apply(beta_diff, c(2, 3), stats::quantile, probs = probs[1],
                 na.rm = TRUE, names = FALSE)
  upper <- apply(beta_diff, c(2, 3), stats::quantile, probs = probs[2],
                 na.rm = TRUE, names = FALSE)

  df <- expand.grid(coef = seq_len(p), diff = seq_len(J - 1),
                    KEEP.OUT.ATTRS = FALSE)
  df$mean <- as.vector(means)
  df$lower <- as.vector(lower)
  df$upper <- as.vector(upper)
  df$coef <- factor(df$coef, levels = seq_len(p))
  df$diff_lab <- paste0("beta[", df$diff + 1, "] - beta[", df$diff, "]")

  if (is.null(xlim)) xlim <- range(df$lower, df$upper, finite = TRUE)
  if (diff(xlim) == 0) xlim <- xlim + c(-1, 1)

  ggplot2::ggplot(df, ggplot2::aes(x = .data$mean, y = .data$coef)) +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed") +
    ggplot2::geom_errorbar(ggplot2::aes(xmin = .data$lower, xmax = .data$upper),
                  orientation = "y", width = 0.2) +
    ggplot2::geom_point(size = 1.6) +
    ggplot2::facet_wrap(
      ~ diff_lab,
      nrow = mfrow[1],
      ncol = mfrow[2],
      labeller = ggplot2::label_parsed
    ) +
    ggplot2::scale_y_discrete(labels = coef_labels) +
    ggplot2::scale_x_continuous(
      breaks = scales::pretty_breaks(n = 6),
      labels = function(x) {
        lab = as.character(x)
        lab[c(1, length(lab))] = ""
        lab
      },
      limits = xlim,
      expand = ggplot2::expansion(mult = c(0.05, 0.05))
    ) +
    ggplot2::labs(x = NULL, y = NULL) +
    ggplot2::theme_bw(base_size = base_size) +
    ggplot2::theme(
      strip.text = ggplot2::element_text(face = "bold", size = base_size + 1),
      axis.text = ggplot2::element_text(family = "Times", size = 15),
      panel.spacing.x = grid::unit(6, "pt"),
      panel.spacing.y = grid::unit(6, "pt"),
      plot.margin = ggplot2::margin(3, 3, 3, 3)
    )
}
