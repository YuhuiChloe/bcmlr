#' Select the number of changepoints with held-out AUC intervals
#'
#' Fit an upper bound on the number of changepoints, count boundaries whose
#' held-out AUC lower confidence limits exceed a threshold, and select the
#' modal count. Refit the entire series only when the selected count is positive.
#'
#' @inheritParams bcmlr
#' @param data A finite numeric matrix with chronological observations in rows
#'   and covariates in columns. Use a one-column matrix for a univariate series.
#' @param thinning Integer held-out observation spacing, at least 2.
#' @param threshold Finite AUC lower-bound cutoff in `[0, 1]`. A boundary is
#'   counted only when its lower bound is strictly greater than this value.
#' @param max_num_cp Positive integer initial number and upper bound.
#' @param sd_beta Positive Gaussian prior variance, forwarded to both stages.
#'
#' @return A list with `num_cp`, the selected count; `max_cps`, the marginal
#'   changepoint modes from the upper-bound fit; `num_cp_dist`, the count for
#'   each retained selection-stage draw; and `probs_num_cp`, a named numeric
#'   vector of empirical proportions for counts `0:max_num_cp`.
#'   The smallest count is selected when modes tie. If `num_cp > 0`, the
#'   list also contains the full-series conditional output from \code{\link{bcmlr}()}.
#'   If `num_cp == 0`, only these four selection summaries are returned;
#'   there is no refit and no artificial boundary is inserted.
#'
#' @details
#' `num_iter` and `num_warmup` apply to each fitting stage. The first stage
#' standardizes training and held-out observations using training statistics.
#' The refit uses full-series standardization and `thinning = 1`.
#' `min_size >= 2 * thinning` is required for the selection stage.
#' This wrapper retains the original strict training-length bound and calls
#' the updated \code{\link{bcmlr}()} sampler, including its numerical corrections.
#' The count proportions summarize the AUC procedure applied to posterior
#' draws; they are not posterior model probabilities from a trans-dimensional
#' sampler. Location and coefficient summaries condition on the selected count.
#'
#' @examples
#' set.seed(2026)
#' X <- matrix(rep(1, 40), ncol = 1)
#' fit <- bcmlr_model_select(X, max_num_cp = 1, min_size = 10,
#'                           num_iter = 12, num_warmup = 6)
#' fit$num_cp
#' fit$probs_num_cp
#' @seealso \code{\link{bcmlr}()}
#' @md
#' @export
bcmlr_model_select <- function(data, init = "even", thinning = 5, prior_beta = "Gaussian", prior_kappa = "default", alpha_f = 0.1, threshold = 0.5, min_size = 30,
                               max_num_cp = 5, sd_beta = 3, num_iter = 5000, num_warmup = 2500, print_progress = FALSE, print_outputs = FALSE){
  if (!is.matrix(data)){
    stop("Please convert your data to a matrix. If you have factor variables, apply a dummification or one-hot encoding.")
  }
  if (!is.numeric(max_num_cp) || length(max_num_cp) != 1L ||
      !is.finite(max_num_cp) || max_num_cp < 1 || max_num_cp != floor(max_num_cp))
    stop("max_num_cp must be a positive integer.", call. = FALSE)
  if (!is.numeric(min_size) || length(min_size) != 1L ||
      !is.finite(min_size) || min_size < 1 || min_size != floor(min_size))
    stop("min_size must be a positive integer.", call. = FALSE)
  if (!is.numeric(threshold) || length(threshold) != 1L ||
      !is.finite(threshold) || threshold < 0 || threshold > 1)
    stop("threshold must be one finite number in [0, 1].", call. = FALSE)
  if (!is.logical(print_outputs) || length(print_outputs) != 1L || is.na(print_outputs))
    stop("print_outputs must be TRUE or FALSE.", call. = FALSE)
  if (max_num_cp > (nrow(data)/min_size -1)){
    stop("The value of max_num_cp exceeded the maximum possible number of changepoints under the specified min_size.")
  }
  # Select the number of changepoints
  old_out = bcmlr(data, init = init, num_CP = max_num_cp, prior_beta = prior_beta, prior_kappa = prior_kappa, alpha_f = alpha_f, thinning = thinning,
              min_size = min_size, sd_beta = sd_beta, num_iter = num_iter, num_warmup = num_warmup,
              print_outputs = FALSE, print_progress = print_progress, model_selection = TRUE)
  if (is.null(old_out$AUC) || any(!is.finite(old_out$AUC)) ||
      !identical(dim(as.matrix(old_out$AUC)),
                 as.integer(c(num_iter - num_warmup, max_num_cp))))
    stop("Unexpected or non-finite AUC output from the selection fit.")
  boolean = (old_out$AUC <= threshold) # TRUE if less than threshold, FALSE otherwise
  rej = rowSums(!as.matrix(boolean)) # Count number of FALSE in each row
  tb = table(factor(rej, levels = seq_len(max_num_cp+1)-1))/length(rej)  # The levels start from 0 to max_num_cp
  num_cp = as.integer(names(tb)[which.max(tb)]) # the posterior mode is taken as the estimated num of CP
  probs = as.numeric(tb) # empirical AUC-selection count proportions from 0 to max_num_cp
  names(probs) = names(tb)
  if (num_cp == 0){
    if (print_outputs) message("There is no change.")
    out = list()
    out$num_cp = num_cp # estimate of the number of CPs
    out$max_cps = old_out$Kappa_mode # the max number of posterior modes from the "old output" before refitting on all the data
    out$num_cp_dist = rej # AUC-selected count at each retained draw
    out$probs_num_cp = probs
  }else{
    # Refit with all the data:
    out = bcmlr(data, init = init, num_CP = num_cp, prior_beta = prior_beta, prior_kappa = prior_kappa, thinning = 1, min_size = min_size, num_iter = num_iter, num_warmup = num_warmup,
                print_outputs = FALSE, sd_beta = sd_beta, print_progress = print_progress, model_selection = FALSE)
    out$num_cp = num_cp # estimate of the number of CPs
    out$max_cps = old_out$Kappa_mode # the max number of posterior modes from the "old output" before refitting on all the data
    out$num_cp_dist = rej # AUC-selected count at each retained draw
    out$probs_num_cp = probs
  }
  if (print_outputs){
    utils::str(out)
  }
  return(out)
}
