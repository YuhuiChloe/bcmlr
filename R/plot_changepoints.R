#' Plot network statistics, changepoint intervals, and optional comparisons
#'
#' @param S Matrix or data frame of network statistics in chronological order.
#' @param dates Date vector aligned with the rows of `S`.
#' @param kappa Posterior changepoint draws: a vector for one changepoint or a
#'   matrix with draws in rows and changepoints in columns.
#' @param other_cps Named list of competitor changepoint-index vectors. Use
#'   `list()` or omit this argument for no comparisons; scalar NA is not accepted.
#' @param kappa_mode Optional changepoint point estimates. NULL computes marginal
#'   posterior modes, choosing the smallest index in a tie.
#' @param stat_names Optional statistic labels in column order.
#' @param method_order Order of the competing methods in the comparison panel.
#' @param ci_level Credible-interval level; default 0.95.
#' @param date_breaks Spacing of date-axis breaks, such as "3 months".
#' @param date_labels Date-label format, such as "%m/%d/%Y".
#' @param title Optional plot title.
#' @param outfile Optional output filename passed to `ggsave()`.
#' @param width,height Plot dimensions in inches when saving.
#' @param dpi Raster resolution when saving.
#' @returns A ggplot object without competitors, or a patchwork object when
#'   a competitor panel is drawn.
#' @details Keeps the paper script's default quantiles and its rounding and
#'   clamping of indices. A changepoint indexes the last network of the preceding
#'   segment. `patchwork` is needed only when competitor estimates are displayed.
#'   Date labels and ticks are restored when no comparison panel is present.
#' @examples
#' S <- cbind(edges = 1:20, triangles = (1:20)^2)
#' dates <- as.Date("2020-01-01") + 7 * (0:19)
#' plot_netstats_cp_cis(S, dates, c(9, 10, 10, 11))
#' @seealso [plot_changepoints()]
#' @export
plot_netstats_cp_cis <- function(S, dates, kappa,
                                 other_cps = list(),
                                 kappa_mode = NULL,
                                 #kappa_joint_mode = NULL, #######
                                 stat_names = NULL,
                                 method_order = names(other_cps),
                                 ci_level = 0.95,
                                 date_breaks = "3 months",
                                 date_labels = "%m/%d/%Y",
                                 title = NULL,
                                 outfile = NULL,
                                 width = 7,
                                 height = 4.8,
                                 dpi = 600) {
  S <- as.data.frame(S)
  Tt <- nrow(S)
  if (length(dates) != Tt) stop("length(dates) must equal nrow(S).")
  if (is.null(stat_names)) stat_names <- colnames(S)
  if (is.null(stat_names)) stat_names <- paste0("Stat ", seq_len(ncol(S)))
  if (length(stat_names) != ncol(S)) stop("stat_names must have one label per column of S.")
  names(S) <- stat_names

  kappa <- as.matrix(kappa)
  if (!is.numeric(kappa) || !length(kappa)) stop("kappa must be numeric posterior changepoint samples.")
  if (!is.null(kappa_mode) && !is.numeric(kappa_mode)) stop("kappa_mode must be numeric when supplied.")
  if (!is.list(other_cps)) stop("other_cps must be a named list of numeric changepoint vectors.")
  if (length(other_cps) && (is.null(names(other_cps)) || any(!nzchar(names(other_cps))))) {
    stop("other_cps must be a named list.")
  }
  if (length(other_cps) && !all(vapply(other_cps, is.numeric, logical(1)))) {
    stop("Each element of other_cps must be a numeric changepoint vector.")
  }

  idx <- function(x) pmax(1L, pmin(Tt, as.integer(round(x))))
  cp_mode <- function(x) {
    x <- idx(x[is.finite(x)])
    if (!length(x)) return(NA_integer_)
    as.integer(names(which.max(table(x))))
  }

  probs <- c((1 - ci_level) / 2, 1 - (1 - ci_level) / 2)
  lbs <- apply(kappa, 2, stats::quantile, probs = probs[1], na.rm = TRUE, names = FALSE)
  ubs <- apply(kappa, 2, stats::quantile, probs = probs[2], na.rm = TRUE, names = FALSE)
  if (is.null(kappa_mode)) kappa_mode <- apply(kappa, 2, cp_mode)

  k_ci <- data.frame(xmin = dates[idx(lbs)], xmax = dates[idx(ubs)])
  k_ci <- k_ci[stats::complete.cases(k_ci), , drop = FALSE]
  k_mode <- data.frame(x = dates[idx(kappa_mode)])
  k_mode <- k_mode[stats::complete.cases(k_mode), , drop = FALSE]
  # Joint posterior mode of the changepoint tuple ####### (two lines)
  #k_joint_mode <- data.frame(x = dates[idx(kappa_joint_mode)])
  #k_joint_mode <- k_joint_mode[stats::complete.cases(k_joint_mode), , drop = FALSE]

  stats_df <- data.frame(
    date = rep(dates, times = ncol(S)),
    series = factor(rep(stat_names, each = Tt), levels = stat_names),
    value = as.vector(as.matrix(S))
  )

  x_scale <- ggplot2::scale_x_date(
    limits = range(dates),
    breaks = seq(min(dates), max(dates), by = date_breaks),
    date_labels = date_labels,
    expand = ggplot2::expansion(mult = c(0.01, 0.02))
  )

  p_stats <- ggplot2::ggplot(stats_df, ggplot2::aes(.data$date, .data$value)) +
    ggplot2::geom_rect(data = k_ci,
                       ggplot2::aes(xmin = .data$xmin, xmax = .data$xmax),
                       inherit.aes = FALSE,
                       ymin = -Inf, ymax = Inf,
                       fill = "steelblue", alpha = 0.18, color = NA) +
    ggplot2::geom_vline(data = k_mode,
                        ggplot2::aes(xintercept = .data$x),
                        linetype = "dashed", linewidth = 0.6, color = "steelblue4") +
    # Joint posterior mode: red dashed lines
    # ggplot2::geom_vline(data = k_joint_mode,
    #                     ggplot2::aes(xintercept = .data$x),
    #                     linetype = "dashed", linewidth = 0.6, color = "firebrick") +
    ggplot2::geom_line(linewidth = 0.4, color = "black") +
    ggplot2::facet_wrap(~ series, ncol = 1, scales = "free_y") +
    x_scale +
    ggplot2::labs(x = NULL, y = NULL) +
    ggplot2::theme_bw(base_family = "serif") +
    ggplot2::theme(
      strip.background = ggplot2::element_rect(fill = "grey95", color = NA),
      # strip.text controls the style of the subplots' titles
      strip.text = ggplot2::element_text(family = "serif", face = "bold", size = 12),
      panel.spacing = grid::unit(6, "pt"),
      axis.text.x = ggplot2::element_blank(),
      axis.ticks.x = ggplot2::element_blank()
    )

  cp_df <- do.call(rbind, lapply(names(other_cps), function(m) {
    cp <- idx(other_cps[[m]][is.finite(other_cps[[m]])])
    cp <- sort(unique(cp))
    if (!length(cp)) return(NULL)
    data.frame(method = m, cp = cp, date = dates[cp])
  }))

  if (!is.data.frame(cp_df) || !nrow(cp_df)) {
    out <- p_stats +
      ggplot2::labs(x = "Date", title = title) +
      ggplot2::theme(
        axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
        axis.ticks.x = ggplot2::element_line()
      )
  } else {
    if (!requireNamespace("patchwork", quietly = TRUE)) {
      stop("Install the suggested package 'patchwork' to show competing methods.")
    }
    cp_df$method <- factor(cp_df$method, levels = method_order)
    p_methods <- ggplot2::ggplot(cp_df, ggplot2::aes(.data$date, .data$method)) +
      ggplot2::geom_point(
        shape = 124,
        size = 4,
        color = "firebrick"
      ) +
      x_scale +
      ggplot2::scale_y_discrete(drop = FALSE) +
      ggplot2::labs(
        x = "Date",
        y = NULL
      ) +
      ggplot2::theme_bw(
        base_size = 13,
        base_family = "serif"
      ) +
      ggplot2::theme(
        # axis.title.x enlarges Date.
        axis.title.x = ggplot2::element_text(
          size = 14,
          margin = ggplot2::margin(t = 10)
        ),
        # axis.text.x enlarges the date labels.
        axis.text.x = ggplot2::element_text(
          size = 11,
          angle = 45,
          hjust = 1
        ),
        # axis.text.y enlarges the competing-method names.
        axis.text.y = ggplot2::element_text(
          size = 12
        )
      )
    # output
    out <- patchwork::wrap_plots(
      p_stats,
      p_methods,
      ncol = 1,
      heights = c(3.2, 1.1)
    ) +
      patchwork::plot_annotation(
        title = title,
        theme = ggplot2::theme(
          plot.title = ggplot2::element_text(
            family = "serif",
            face = "bold",
            size = 15,
            hjust = 0.5
          )
        )
      )
  }

  if (!is.null(outfile)) ggplot2::ggsave(outfile, out, width = width, height = height, dpi = dpi)
  return(out)
}


#' Plot changepoints using the dissertation's argument order
#'
#' A thin wrapper around [plot_netstats_cp_cis()] that accepts the posterior
#' draws as its first argument. All computations and plotting settings are
#' delegated to that function.
#' @inheritParams plot_netstats_cp_cis
#' @param ... Other arguments passed to [plot_netstats_cp_cis()], including
#'   `other_cps`, `kappa_mode`, and `stat_names`.
#' @returns The ggplot or patchwork object returned by [plot_netstats_cp_cis()].
#' @examples
#' S <- cbind(edges = 1:20, triangles = (1:20)^2)
#' dates <- as.Date("2020-01-01") + 7 * (0:19)
#' plot_changepoints(c(9, 10, 10, 11), S = S, dates = dates)
#' @export
plot_changepoints <- function(kappa, S, dates, ...) {
  plot_netstats_cp_cis(S = S, dates = dates, kappa = kappa, ...)
}
