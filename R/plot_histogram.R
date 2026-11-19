#' Plot the posterior distribution of changepoint count
#'
#' @param num_vec Posterior draws of the number of changepoints, usually
#'   `fit$num_cp_dist`. As in the paper script, values are converted to integers
#'   and missing values are removed. Use nonnegative integer counts.
#' @param max_num Largest displayed count; NULL uses the largest retained count.
#' @param title,subtitle Plot title and optional subtitle.
#' @param font_size Size of the probability labels, passed to `geom_text()`.
#' @param xlab,ylab Axis labels.
#' @param y_max Optional upper y-axis limit; NULL chooses it automatically.
#' @param highlight_mode Highlight the smallest modal count in a tie?
#' @returns A ggplot object with empirical posterior probability bars.
#' @examples
#' plot_histogram(c(2, 3, 3, 3, 4), max_num = 5)
#' @export
plot_histogram <- function(num_vec,
                            max_num = NULL,
                            title = "Posterior distribution of changepoint count",
                            subtitle = NULL,
                            font_size = 4,
                            xlab = "Number of changepoints",
                            ylab = "Posterior probability",
                            y_max = NULL,
                            highlight_mode = TRUE) {

  # Clean input
  num_vec <- as.integer(num_vec)
  num_vec <- num_vec[!is.na(num_vec)]

  if (length(num_vec) == 0) {
    stop("num_vec contains no non-missing values.")
  }

  if (is.null(max_num)) {
    max_num <- max(num_vec)
  }

  if (max_num < max(num_vec)) {
    stop("max_num must be at least max(num_vec).")
  }

  # Posterior probability mass table
  tab <- as.data.frame(table(factor(num_vec, levels = 0:max_num)))
  names(tab) <- c("num_cp", "count")

  tab$num_cp <- as.integer(as.character(tab$num_cp))
  tab$prob <- tab$count / sum(tab$count)
  tab$label <- sprintf("%.3f", tab$prob)
  tab$y_text <- ifelse(tab$prob == 0, 0.01, tab$prob)

  mode_cp <- tab$num_cp[which.max(tab$prob)]
  tab$is_mode <- tab$num_cp == mode_cp

  # if (is.null(subtitle)) {
  #   subtitle <- paste0(
  #     "Posterior mode = ", mode_cp,
  #     "; Pr(K = ", mode_cp, ") = ",
  #     sprintf("%.3f", max(tab$prob))
  #   )
  # }

  if (is.null(y_max)) {
    y_max <- min(1, max(tab$prob) + 0.12)
    y_max <- max(y_max, 0.25)
  }

  ggplot2::ggplot(tab, ggplot2::aes(x = .data$num_cp, y = .data$prob)) +
    ggplot2::geom_col(
      ggplot2::aes(fill = if (highlight_mode) .data$is_mode else FALSE),
      width = 0.72,
      color = "grey30",
      linewidth = 0.25,
      show.legend = FALSE
    ) +
    ggplot2::geom_text(
      ggplot2::aes(y = .data$y_text,label = .data$label),
      vjust = -0.35,
      size = font_size
    ) +
    ggplot2::scale_fill_manual(
      values = c("FALSE" = "grey82", "TRUE" = "#4C6A92")
    ) +
    ggplot2::scale_x_continuous(
      breaks = 0:max_num,
      expand = ggplot2::expansion(add = 0.35)
    ) +
    ggplot2::scale_y_continuous(
      breaks = scales::pretty_breaks(n = 6),
      labels = scales::label_number(accuracy = 0.1),
      limits = c(0, y_max),
      expand = ggplot2::expansion(mult = c(0, 0.06))
    ) +
    ggplot2::labs(
      title = title,
      subtitle = subtitle,
      x = xlab,
      y = ylab
    ) +
    ggplot2::theme_classic(base_size = 12, base_family = "serif") +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", size = 16, hjust = 0.5),
      plot.subtitle = ggplot2::element_text(hjust = 0.5, size = 15),
      axis.text = ggplot2::element_text(color = "black"),
      axis.title.x = ggplot2::element_text(margin = ggplot2::margin(t = 15)),
      axis.title.y = ggplot2::element_text(margin = ggplot2::margin(r = 15))
    )
}
