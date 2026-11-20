#' Plot posterior intervals for adjacent coefficient differences
#'
#' @param Beta Numeric array of posterior draws by network statistics by
#'   chronological segments, as returned in `fit$Beta`. A draws-by-covariates
#'   matrix from a one-changepoint bcmlr fit is also accepted; it is the later
#'   segment's coefficient relative to the zero first segment.
#' @param mfrow Optional panel layout, c(rows, columns). NULL chooses a layout
#'   with at most three columns and enough panels for every changepoint.
#' @param coef_labels Optional labels in the same order as the columns of `S`.
#' @param ci_level Credible-interval level; default 0.95.
#' @param xlim Optional finite increasing x-axis limits; zooms the display
#'   without dropping posterior summaries.
#' @param base_size Base text size for `theme_bw()`.
#' @returns A ggplot object. Points show posterior means of draw-wise
#'   next-minus-previous coefficient differences; bars show equal-tailed intervals.
#' @details Differences are calculated within each draw. The one-changepoint
#'   matrix already represents the later-minus-earlier contrast. Serif fonts
#'   and automatic panel layout are used consistently across interpretation plots.
#' @examples
#' Beta <- array(0, c(5, 2, 2))
#' Beta[, 1, 2] <- seq(0.1, 0.5, length.out = 5)
#' plot_beta_diff_CIs(Beta, mfrow = c(1, 1), coef_labels = c("Edges", "Triangles"))
#' @export
plot_beta_diff_CIs <- function(Beta, mfrow = NULL, coef_labels = NULL,
                               ci_level = 0.95, xlim = NULL,
                               base_size = 12) {
  # The single-changepoint sampler uses a zero first reference class.
  if (!is.numeric(Beta) || is.complex(Beta) ||
      !length(dim(Beta)) %in% c(2L, 3L) || any(dim(Beta) < 1L) ||
      any(!is.finite(Beta))) stop("Beta must be a finite numeric matrix or array.")
  if (is.matrix(Beta)) {
    binary <- array(0, c(nrow(Beta), ncol(Beta), 2L))
    binary[, , 2L] <- Beta
    Beta <- binary
  }
  if (!is.numeric(ci_level) || length(ci_level) != 1L || !is.finite(ci_level) ||
      ci_level <= 0 || ci_level >= 1) stop("ci_level must lie in (0, 1).")

  p <- dim(Beta)[2]
  J <- dim(Beta)[3]
  if (J < 2) stop("Beta must have at least two classes.")
  if (is.null(mfrow)) {
    nc <- min(3L, J - 1L)
    mfrow <- c(ceiling((J - 1L) / nc), nc)
  }
  if (!is.numeric(mfrow) || length(mfrow) != 2L || any(!is.finite(mfrow)) ||
      any(mfrow < 1 | mfrow != floor(mfrow)) || prod(mfrow) < J - 1L) {
    stop("mfrow must provide enough panels for every changepoint.")
  }
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
  labels <- paste0("beta[", seq_len(J - 1L) + 1L, "] - beta[", seq_len(J - 1L), "]")
  df$diff_lab <- factor(labels[df$diff], levels = labels)

  if (is.null(xlim)) xlim <- range(df$lower, df$upper, finite = TRUE)
  if (!is.numeric(xlim) || length(xlim) != 2L || any(!is.finite(xlim)) ||
      xlim[1] > xlim[2]) stop("xlim must contain two finite ordered limits.")
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
      expand = ggplot2::expansion(mult = c(0.05, 0.05))
    ) +
    ggplot2::coord_cartesian(xlim = xlim) +
    ggplot2::labs(x = NULL, y = NULL) +
    ggplot2::theme_bw(base_size = base_size, base_family = "serif") +
    ggplot2::theme(
      strip.text = ggplot2::element_text(face = "bold", size = base_size + 1),
      axis.text = ggplot2::element_text(family = "serif", size = base_size - 1),
      panel.spacing.x = grid::unit(6, "pt"),
      panel.spacing.y = grid::unit(6, "pt"),
      plot.margin = ggplot2::margin(3, 3, 3, 3)
    )
}
