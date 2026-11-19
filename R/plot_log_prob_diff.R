#' Plot posterior log-probability differences between adjacent segments
#'
#' Summarize draw-wise log-probability differences over the two segments
#' neighboring each changepoint, as in Chapter 4 of the dissertation.
#'
#' @param P Posterior class-probability draws from a full-series conditional
#'   bcmlr fit. For multiple changepoints, a numeric array with dimensions
#'   draws by observations by chronological segments. For one changepoint,
#'   a numeric matrix with dimensions draws by observations, containing the
#'   probability of the second segment; a two-segment array is also accepted.
#' @param kappa_mode Integer changepoint point estimates, usually
#'   \code{fit$Kappa_mode}. They index the last observation of each preceding
#'   segment and must form a strictly increasing partition of the rows in P.
#' @param mfrow Optional panel layout, c(rows, columns). NULL chooses a layout
#'   with at most three columns and enough panels for every changepoint.
#' @param dates Optional strictly increasing Date vector, one per observation.
#'   NULL displays observation indices.
#' @param ci_level Pointwise equal-tailed credible-interval level; default 0.95.
#' @param base_size Base text size for \code{theme_bw()}.
#' @param date_breaks,date_labels Date-axis spacing and label format.
#' @param title Optional plot title.
#' @param outfile Optional output filename passed to \code{ggsave()}.
#' @param width,height Plot dimensions in inches when saving.
#' @param dpi Raster resolution when saving.
#' @param eps Optional probability floor strictly between 0 and 0.5. NULL
#'   computes exact log differences and rejects zero probabilities (also ones
#'   in the single-changepoint matrix). An explicit floor allows approximate
#'   summaries of saturated draws and emits a warning whenever values change.
#'   For a matrix, the second-segment probability is clipped to
#'   \code{[eps, 1 - eps]}; for an array, each class probability is floored.
#' @returns A ggplot object. Its \code{data} contains the changepoint number
#'   \code{l}, observation index \code{idx}, plotted \code{time}, posterior
#'   \code{mean}, \code{lower} and \code{upper} interval limits, changepoint
#'   index \code{cp}, and ordered facet label \code{lab}.
#' @details For each retained draw, this function computes
#'   \code{log(q_next) - log(q_current)} before taking means and quantiles.
#'   It does not compute a log ratio of posterior mean probabilities. Taking
#'   the two logs separately avoids overflow in an intermediate ratio.
#'
#'   With one changepoint, the contrast is \code{log(P) - log1p(-P)}, and the
#'   panel covers the entire series. More generally, panel l covers indices
#'   \code{kappa[l - 1] + 1} through \code{kappa[l + 1]}, with boundary values
#'   zero and the series length. Bands are pointwise credible intervals.
#'
#'   Probabilities in each draw and observation of an array must sum to one
#'   within 1e-8. Nonfinite draws are rejected, not silently discarded.
#'   If marginal modes do not form a valid partition, supply a valid joint
#'   changepoint draw explicitly; estimates are never sorted or clamped.
#'   Use the full-data refit returned by \code{bcmlr_model_select()}, or a
#'   conditional \code{bcmlr()} fit with \code{thinning = 1}. A thinned fit's
#'   probability columns and original-series changepoint indices are not aligned.
#' @examples
#' P <- matrix(c(0.2, 0.3, 0.4, 0.6, 0.7, 0.8), nrow = 2)
#' plot_log_prob_diff(P, kappa_mode = 1)
#' @seealso \code{\link{bcmlr}}, \code{\link{bcmlr_model_select}}
#' @importFrom rlang .data
#' @export
plot_log_prob_diff <- function(P, kappa_mode, mfrow = NULL, dates = NULL,
                               ci_level = 0.95, base_size = 12,
                               date_breaks = "6 months",
                               date_labels = "%m/%d/%Y", title = NULL,
                               outfile = NULL, width = 10, height = 3.8,
                               dpi = 600, eps = NULL) {
  # Validate the two probability formats returned by bcmlr.
  dims <- dim(P)
  if (!is.numeric(P) || is.complex(P) || !length(dims) %in% c(2L, 3L) ||
      any(dims < 1L) || dims[2L] < 2L || any(!is.finite(P)) ||
      any(P < 0 | P > 1)) {
    stop("P must be a finite probability matrix or array with at least two observations.",
         call. = FALSE)
  }
  binary <- length(dims) == 2L
  n_draws <- dims[1L]
  n_time <- dims[2L]
  num_cp <- if (binary) 1L else dims[3L] - 1L
  if (num_cp < 1L) stop("There are no changepoints to interpret.", call. = FALSE)
  if (!binary && any(abs(apply(P, c(1, 2), sum) - 1) > 1e-8)) {
    stop("Class probabilities must sum to one for every draw and observation.",
         call. = FALSE)
  }
  if (!is.numeric(kappa_mode) || is.complex(kappa_mode) ||
      length(kappa_mode) != num_cp || any(!is.finite(kappa_mode)) ||
      any(kappa_mode != floor(kappa_mode)) ||
      any(kappa_mode < 1 | kappa_mode >= n_time) ||
      any(diff(kappa_mode) <= 0)) {
    stop("kappa_mode must be a strictly increasing valid partition aligned with P; supply a valid joint draw if needed.",
         call. = FALSE)
  }
  if (!is.numeric(ci_level) || is.complex(ci_level) || length(ci_level) != 1L ||
      !is.finite(ci_level) || ci_level <= 0 || ci_level >= 1) {
    stop("ci_level must lie strictly between zero and one.", call. = FALSE)
  }
  if (!is.numeric(base_size) || is.complex(base_size) || length(base_size) != 1L ||
      !is.finite(base_size) || base_size <= 0) {
    stop("base_size must be a positive finite number.", call. = FALSE)
  }
  if (!is.null(dates) && (!inherits(dates, "Date") || length(dates) != n_time ||
      any(!is.finite(as.numeric(dates))) || any(diff(dates) <= 0))) {
    stop("dates must be a strictly increasing Date vector aligned with P.", call. = FALSE)
  }
  if (is.null(mfrow)) {
    nc <- min(3L, num_cp)
    mfrow <- c(ceiling(num_cp / nc), nc)
  }
  if (!is.numeric(mfrow) || is.complex(mfrow) || length(mfrow) != 2L ||
      any(!is.finite(mfrow)) || any(mfrow < 1 | mfrow != floor(mfrow)) ||
      prod(mfrow) < num_cp) {
    stop("mfrow must provide enough panels for every changepoint.", call. = FALSE)
  }

  # Clipping is opt-in because it changes the posterior summaries.
  if (!is.null(eps)) {
    if (!is.numeric(eps) || is.complex(eps) || length(eps) != 1L ||
        !is.finite(eps) || eps <= 0 || eps >= 0.5 || (binary && 1 - eps == 1)) {
      stop("eps must lie in (0, 0.5); for a matrix, 1 - eps must be smaller than 1.",
           call. = FALSE)
    }
    changed <- P < eps | (binary & P > 1 - eps)
    if (any(changed)) {
      warning("Probabilities were clipped using eps; the plotted log differences are approximate.",
              call. = FALSE)
      P[P < eps] <- eps
      if (binary) P[P > 1 - eps] <- 1 - eps
    }
  }
  if (any(P <= 0) || (binary && any(P >= 1))) {
    stop("Zero class probabilities give infinite log differences. Supply an explicit eps only if approximate clipping is intended.",
         call. = FALSE)
  }

  # Compute each contrast draw by draw on its two neighboring segments.
  bounds <- c(0L, kappa_mode, n_time)
  time <- if (is.null(dates)) seq_len(n_time) else dates
  probs <- c((1 - ci_level) / 2, (1 + ci_level) / 2)
  labels <- paste0("bold(x)[i]^T * (bold(beta)[", seq_len(num_cp) + 1L,
                   "] - bold(beta)[", seq_len(num_cp), "])")
  panels <- lapply(seq_len(num_cp), function(l) {
    idx <- seq.int(bounds[l] + 1L, bounds[l + 2L])
    if (binary) {
      q_next <- P[, idx, drop = FALSE]
      differences <- log(q_next) - log1p(-q_next)
    } else {
      q_next <- matrix(P[, idx, l + 1L, drop = FALSE], nrow = n_draws)
      q_current <- matrix(P[, idx, l, drop = FALSE], nrow = n_draws)
      differences <- log(q_next) - log(q_current)
    }
    intervals <- vapply(seq_along(idx), function(i) {
      stats::quantile(differences[, i], probs = probs, names = FALSE)
    }, numeric(2))
    data.frame(l = l, idx = idx, time = time[idx],
               mean = colMeans(differences), lower = intervals[1, ],
               upper = intervals[2, ], cp = kappa_mode[l], lab = labels[l])
  })
  df <- do.call(rbind, panels)
  df$lab <- factor(df$lab, levels = labels)
  cp_df <- data.frame(time = time[kappa_mode], lab = factor(labels, levels = labels))

  # Match the ggplot, serif font, and faceting style used in bcmlrNET.
  p <- ggplot2::ggplot(df, ggplot2::aes(x = .data$time, y = .data$mean)) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed") +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = .data$lower, ymax = .data$upper),
                         fill = "blue", alpha = 0.25) +
    ggplot2::geom_line(linewidth = 0.6) +
    ggplot2::geom_point(size = 1.4) +
    ggplot2::geom_vline(data = cp_df, ggplot2::aes(xintercept = .data$time),
                        color = "blue", linetype = "dashed", linewidth = 0.8) +
    ggplot2::facet_wrap(~ lab, nrow = mfrow[1], ncol = mfrow[2],
                         labeller = ggplot2::label_parsed, scales = "free_x") +
    ggplot2::labs(x = if (is.null(dates)) "Observation index" else "Date",
                   y = NULL, title = title) +
    ggplot2::theme_bw(base_size = base_size, base_family = "serif") +
    ggplot2::theme(strip.text = ggplot2::element_text(face = "bold", size = base_size + 1),
                   axis.text = ggplot2::element_text(size = base_size - 1),
                   panel.spacing.x = grid::unit(6, "pt"),
                   panel.spacing.y = grid::unit(6, "pt"),
                   plot.margin = ggplot2::margin(6, 10, 6, 6))
  if (is.null(dates)) {
    p <- p + ggplot2::scale_x_continuous(breaks = scales::pretty_breaks(n = 6),
                                         expand = ggplot2::expansion(mult = c(0.02, 0.02)))
  } else {
    p <- p + ggplot2::scale_x_date(date_breaks = date_breaks, date_labels = date_labels,
                                   expand = ggplot2::expansion(mult = c(0.02, 0.02)))
  }
  if (!is.null(outfile)) ggplot2::ggsave(outfile, p, width = width, height = height, dpi = dpi)
  return(p)
}
