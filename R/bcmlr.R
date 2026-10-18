#' Generalized Bayesian changepoint analysis for a fixed count
#'
#' Fit the multinomial-log-loss changepoint model using Polya-Gamma Gibbs
#' updates and Gaussian or horseshoe coefficient priors. This implementation
#' retains separate single- and multiple-changepoint branches and their
#' historical output formats.
#'
#' @param data A finite numeric matrix with observations in chronological rows
#'   and covariates in columns. A numeric vector is accepted as a univariate
#'   series. No intercept is added. Constant columns are centered to zero.
#' @param num_CP A positive integer number of changepoints. To allow selection
#'   of zero changepoints, use \code{\link{bcmlr_model_select}()}.
#' @param init One of `"even"`, `"ecp"`, `"kcp"`, or `"multirank"`.
#'   A single changepoint is initialized at the training-series midpoint.
#'   For multiple changepoints, `"even"` uses approximately equal segments.
#'   `"ecp"` requires ecp; the other two external methods require reticulate
#'   and Python ruptures. External estimates must satisfy the fitting support.
#' @param prior_beta Either `"Gaussian"` or `"horseshoe"`.
#' @param prior_kappa `"default"` for a prior proportional to the product of
#'   `segment_length^(-segment_length)`, or a function returning one finite,
#'   nonnegative joint prior mass for the ordered changepoint vector. The
#'   function receives training-series indices when observations are held out.
#' @param alpha_f AUC confidence-interval error probability in `(0, 1)`.
#'   The default 0.1 gives 90 percent intervals. Used for model selection only.
#' @param min_size Positive integer strict lower bound on segment lengths.
#'   With held-out spacing `k > 1`, the bound becomes
#'   `round(min_size * (1 - 1/k))` on the training series. Every training
#'   segment must have length strictly greater than that bound.
#' @param sd_beta Positive Gaussian prior variance, despite the historical
#'   name. The prior covariance is `sd_beta * I`. Ignored by the horseshoe fit.
#' @param tempering Positive power on the log-loss contribution, with the
#'   priors left untempered. Use 1 for ordinary inference. This argument does
#'   not implement parallel tempering or chain swaps. Values below 1 are not
#'   supported by the selected Polya-Gamma sampler and are rejected.
#' @param thinning Positive integer held-out observation spacing, not thinning
#'   of MCMC draws. At 1 all observations are used for fitting; at `k > 1`,
#'   every kth original observation is held out.
#' @param num_iter Total number of MCMC iterations, a positive integer.
#' @param num_warmup Number of initial iterations discarded, a nonnegative
#'   integer strictly smaller than `num_iter`.
#' @param print_outputs Whether to print the returned object's structure.
#' @param print_progress Whether to display a progress bar.
#' @param model_selection Whether to calculate held-out AUC lower confidence
#'   limits after warmup. Requires `thinning >= 2`, `tempering = 1`, and
#'   `min_size >= 2 * thinning`. The actual held-out class sizes are checked.
#'
#' @return A list. Let `S = num_iter - num_warmup`, `p = ncol(data)`,
#'   `N` be the number of training observations, and `J = num_CP + 1`.
#'   For a conditional fit (`model_selection = FALSE`):
#'   * `Kappa`: an S-vector for one changepoint, or an `S x num_CP` matrix.
#'     Entries are last-observation indices in the original series.
#'   * `Kappa_mode`: marginal posterior modes, with the smallest mode in a tie.
#'   * `Beta`: an `S x p` matrix for one changepoint; the first segment is
#'     the zero reference class, so coefficients describe the later segment.
#'     For multiple changepoints, an `S x p x J` array with the last class zero.
#'   * `P`: an `S x N` matrix of later-class probabilities for one changepoint,
#'     or an `S x N x J` array for multiple changepoints. When `thinning > 1`,
#'     these probabilities cover training rows only, unlike `Kappa` indices.
#'   * `Beta_mean`, `P_mean`: corresponding posterior means.
#'
#'   For a selection-stage fit (`model_selection = TRUE`), returns `Kappa`,
#'   `Kappa_mode`, and `AUC`; `AUC` is an S-vector for one changepoint or an
#'   `S x num_CP` matrix of held-out AUC lower confidence limits. Coefficient
#'   and probability draws are not returned in this mode.
#'
#' @details
#' Centering and standardization use training means and sample standard
#' deviations, also applied to held-out observations. The original compressed
#' training-series support and prior-length convention are retained. A custom
#' prior must be evaluable without relying on undeclared global objects.
#'
#' The horseshoe hierarchy uses one global scale shared by all nonreference
#' coefficients; its global conditional therefore uses `p * num_CP` coefficients.
#' For one changepoint this reduces to p coefficients. Local scales are separate.
#'
#' This revision repairs detached default-prior terms, horseshoe precision and
#' auxiliary updates, and ROC orientation. These corrections can change results
#' from earlier scripts. They are not an assertion of exact reproduction of
#' historical analyses. Fixed ROC direction treats the later segment as positive.
#' DeLong intervals for perfectly separated small samples can be degenerate;
#' pROC warnings are retained. Inference is conditional on `num_CP`.
#'
#' @examples
#' set.seed(2026)
#' X <- matrix(c(rnorm(20), rnorm(20, 2)), ncol = 1)
#' fit <- bcmlr(X, num_CP = 1, min_size = 5,
#'              num_iter = 20, num_warmup = 10)
#' fit$Kappa_mode
#' dim(fit$Beta)
#' @seealso \code{\link{bcmlr_model_select}()}
#' @md
#' @export
bcmlr <- function(data, num_CP, init = "even", prior_beta = "Gaussian", prior_kappa = "default", alpha_f = 0.1,
                  min_size = 30, sd_beta = 3, tempering = 1, thinning = 1,
                  num_iter = 5000, num_warmup = 2500, print_outputs = FALSE,
                  print_progress = FALSE, model_selection = FALSE){
  # Validate scalars before comparing them or constructing any sequences.
  integer_scalar <- function(x, name, lower = 1) {
    if (!is.numeric(x) || length(x) != 1L || !is.finite(x) ||
        x < lower || x != floor(x) || x > .Machine$integer.max) {
      stop(name, " must be a finite integer >= ", lower, ".", call. = FALSE)
    }
  }
  flag <- function(x, name) {
    if (!is.logical(x) || length(x) != 1L || is.na(x))
      stop(name, " must be TRUE or FALSE.", call. = FALSE)
  }
  integer_scalar(num_CP, "num_CP")
  integer_scalar(min_size, "min_size")
  integer_scalar(thinning, "thinning")
  integer_scalar(num_iter, "num_iter")
  integer_scalar(num_warmup, "num_warmup", 0)
  if (num_warmup >= num_iter)
    stop("num_warmup must be smaller than num_iter.", call. = FALSE)
  flag(print_outputs, "print_outputs")
  flag(print_progress, "print_progress")
  flag(model_selection, "model_selection")
  if (!is.numeric(sd_beta) || length(sd_beta) != 1L ||
      !is.finite(sd_beta) || sd_beta <= 0)
    stop("sd_beta must be one positive finite variance.", call. = FALSE)
  if (!is.numeric(alpha_f) || length(alpha_f) != 1L ||
      !is.finite(alpha_f) || alpha_f <= 0 || alpha_f >= 1)
    stop("alpha_f must lie strictly between zero and one.", call. = FALSE)
  if (!is.numeric(tempering) || length(tempering) != 1L ||
      !is.finite(tempering) || tempering < 1)
    stop("tempering must be finite and >= 1 for this Polya-Gamma sampler.",
         call. = FALSE)
  if (!is.character(prior_beta) || length(prior_beta) != 1L ||
      is.na(prior_beta) || !prior_beta %in% c("Gaussian", "horseshoe"))
    stop("prior_beta must be 'Gaussian' or 'horseshoe'.", call. = FALSE)
  if (!identical(prior_kappa, "default") && !is.function(prior_kappa))
    stop("prior_kappa must be 'default' or a joint-prior function.", call. = FALSE)
  if (!is.character(init) || length(init) != 1L || is.na(init) ||
      !init %in% c("even", "kcp", "ecp", "multirank"))
    stop("init must be 'even', 'kcp', 'ecp', or 'multirank'.", call. = FALSE)
  data <- as.matrix(data)
  if (!is.numeric(data) || nrow(data) < 2L || ncol(data) < 1L ||
      any(!is.finite(data)))
    stop("data must contain finite numeric values, at least two rows, and one column.",
         call. = FALSE)
  storage.mode(data) <- "double"
  if (thinning > 1 && thinning >= min_size)
    stop("min_size must exceed thinning when observations are held out.", call. = FALSE)
  if (model_selection && (thinning < 2 || min_size < 2 * thinning || tempering != 1))
    stop("Selection requires thinning >= 2, min_size >= 2 * thinning, and tempering = 1.",
         call. = FALSE)
  if (thinning != 1) min_size <- round(min_size * (1 - 1/thinning))
  holdout_idx <- integer()
  #### Center & standardize the data ####
  # A function to center & standardize a matrix x
  center_standardize <- function(x, m, sdv){
    ans <- sweep(sweep(x, 2L, m, "-"), 2L, sdv, "/")
    if (any(!is.finite(ans)))
      stop("Standardization produced non-finite values; rescale the data.")
    ans
  }
  X_all = as.matrix(data)
  n = nrow(X_all) # the sample size of the input data
  #### Hold out samples ####
  if (thinning != 1){
    all_idx = seq_len(n) # all indices from the input data
    holdout_idx = all_idx[all_idx%%thinning == 0] # the indices of the held-out samples
    train_idx = all_idx[-holdout_idx] # the indices for fitting the model (the trainng set indices)
    X_auc = X_all[holdout_idx, , drop = FALSE] # rows of X selected to calculate AUC
    X = X_all[-holdout_idx, , drop = FALSE] # training observations
    if (nrow(X) < 2L) stop("At least two training observations are required.")
    if (ncol(X) == 1){
      m = mean(X, na.rm = TRUE)
    }else{
      m = colMeans(X, na.rm = TRUE) # column means in the training set
    }
    sdv = apply(X, 2, FUN = stats::sd, na.rm = TRUE) # standard deviation of the training set
    sdv[sdv == 0] = 1 # Avoid divide-by-zero in standardization: sdv = 0 could occur if the input data has a column of the same constant.
    X = center_standardize(X, m, sdv) # center & standardize the subseries used to fit the model with its column means
    X_auc = center_standardize(X_auc, m, sdv) # center & standardize the held-out samples with the same column means
    # Update the original data set with the centered-standardized X and X_auc
    if (length(holdout_idx) > 0) {
      X_all[holdout_idx, ] = X_auc
    }
    if (length(train_idx) > 0) {
      X_all[train_idx,   ] = X
    }
    # Sanity checks
    stopifnot(
      nrow(X_auc) == length(holdout_idx),
      nrow(X)     == length(train_idx)
    )
    if (length(holdout_idx) > 0) {
      stopifnot(isTRUE(all.equal(X_all[holdout_idx, , drop = FALSE], X_auc)))
    }
    if (length(train_idx) > 0) {
      stopifnot(isTRUE(all.equal(X_all[train_idx,   , drop = FALSE], X)))
    }
  }else if (thinning == 1){ # if thinning = 1
    all_idx = seq_len(n)
    # holdout_idx = all_idx
    train_idx = all_idx
    if (dim(X_all)[2] == 1){
      m = mean(X_all, na.rm = TRUE)
    }else{
      m = colMeans(X_all, na.rm = TRUE) # column means in the training set
    }
    sdv = apply(X_all, 2, FUN = stats::sd)
    sdv[sdv == 0] = 1 # Avoid divide-by-zero in standardization: sdv = 0 could occur if the input data has a column of the same constant.
    X = center_standardize(X_all, m, sdv)
    X_all = X
    # No need for X_auc since we are not computing AUCs for model selection
  }
  N = dim(X)[1] # Sample size of training set
  p = dim(X)[2] # Dimension of the training set
  if ((num_CP + 1) * (min_size + 1) > N)
    stop("Too many changepoints for the strict training-series min_size bound.")
  if (num_CP > 1 && init %in% c("kcp", "multirank")) {
    if (!requireNamespace("reticulate", quietly = TRUE))
      stop("Install reticulate and Python ruptures for this initialization.")
    rpt <- reticulate::import("ruptures", delay_load = FALSE)
  }
  if (num_CP > 1 && init == "ecp" && !requireNamespace("ecp", quietly = TRUE))
    stop("Install ecp for init = 'ecp'.")
  #### Utility Functions ####
  # https://arxiv.org/pdf/2205.04997
  # https://gregorygundersen.com/blog/2020/02/09/log-sum-exp/
  # A function that mitigates computational issues in massive division or multiplication of extreme values.
  logsumexp <- function(x) {
    a <- max(x)
    if (identical(a, -Inf)) return(-Inf)
    if (!is.finite(a) || anyNA(x)) stop("Non-finite log weights or predictors.")
    a + log(sum(exp(x - a)))
  }
  log_probabilities <- function(eta) {
    if (any(!is.finite(eta))) stop("Non-finite linear predictor; rescale inputs.")
    eta - apply(eta, 1L, logsumexp)
  }
  positive_finite <- function(x, name) {
    if (any(!is.finite(x)) || any(x <= 0))
      stop(name, " must remain positive and finite.")
  }
  solve_chol <- function(U, b) backsolve(U, forwardsolve(t(U), b))
  candidate_range <- function(left, right) {
    lower <- left + 1 + min_size
    upper <- right - 1 - min_size
    if (lower > upper) stop("Empty changepoint support; check initialization and min_size.")
    seq.int(lower, upper)
  }
  normalize_weights <- function(w) {
    if (anyNA(w) || any(w == Inf) || !any(is.finite(w)))
      stop("No finite positive conditional mass; check the data and changepoint prior.")
    exp(w - logsumexp(w))
  }
  joint_prior_mass <- function(kappa) {
    value <- prior_kappa(kappa)
    if (!is.numeric(value) || length(value) != 1L ||
        !is.finite(value) || value < 0)
      stop("prior_kappa must return one finite nonnegative joint mass.")
    value
  }
  auc_interval <- function(y, score) {
    y <- as.integer(y)
    score <- as.numeric(score)
    if (length(y) != length(score) || any(!is.finite(score)) ||
        anyNA(y) || !all(y %in% c(0L, 1L)) ||
        sum(y == 0L) < 2L || sum(y == 1L) < 2L)
      stop("AUC intervals require two finite held-out observations per class.")
    if (length(unique(score)) == 1L) return(rep(0.5, 3L))
    roc <- pROC::roc(y, score, levels = c(0, 1), direction = "<", quiet = TRUE)
    ans <- as.numeric(pROC::ci.auc(roc, conf.level = 1 - alpha_f,
                                   method = "delong"))
    if (length(ans) != 3L || any(!is.finite(ans))) stop("Non-finite AUC interval.")
    ans
  }
  valid_initialization <- function(kappa) {
    length(kappa) == num_CP && all(is.finite(kappa)) &&
      all(kappa == floor(kappa)) && all(diff(c(0, kappa, N)) > min_size)
  }
  #### Utility functions ####
  # 1-1 Correspondence
  # A function that converts a change point vector to a response matrix with the change point structure.
  kappa_to_Y <- function(kappa, N, L){
    class <- rep(1, kappa[1])
    for(j in 2:L){
      class <- c(class, rep(j, kappa[j] - kappa[j-1]))
    }
    class <- c(class, rep(L+1, N - kappa[L]))
    class <- as.factor(class)
    Y <- stats::model.matrix(~ -1 + class, data=class)
    return(Y)
  }
  # A function that updates \beta, \kappa, and  \omega per iteration  #
  GS_update <- function(temper, beta, kappa, omega, tau_sq, xi, nu, lambda_sq, X_all, X, train_idx, holdout_idx, n, N, p, L, J, V0inv, V0, model_selection){
    ci <- NULL # No AUC is calculated during warmup.
    # 0. Update Y, delta
    Y = kappa_to_Y(kappa, N, L)
    delta = temper*(Y-1/2)
    # 1. Update HS prior (on beta) parameters & beta & omega
    if (prior_beta == "horseshoe"){
      # 1.1 Update HS priors
      xi = LearnBayes::rigamma(n =1, a = 1, b = 1+1/tau_sq)
      xi_inv = 1/xi
      # One global tau is shared across all p * (J - 1) coefficients.
      tau_sq <- LearnBayes::rigamma(n = 1, a = 0.5 * (p * (J - 1) + 1),
                         b = xi_inv + 0.5 * sum(beta[, seq_len(J - 1), drop = FALSE]^2 / lambda_sq))
      for (j in 1:(J-1)){
        for (d in 1:p){
          nu[d,j] <- LearnBayes::rigamma(n=1, a = 1, b = 1+1/lambda_sq[d,j])
          lambda_sq[d,j] <- LearnBayes::rigamma(n=1, a = 1, b = 1/nu[d,j] + 0.5/ tau_sq*beta[d,j]^2)
        }
      }
      # 1.2 Update beta and omega
      for(j in 1:(J-1)){
        positive_finite(lambda_sq[,j] * tau_sq, "Horseshoe variance")
        Sigma0_j = diag(1/(lambda_sq[,j]*tau_sq), nrow = p, ncol = p)
        # Polson et al. (2013) describes or implies the existence of matrices
        # C and eta. We compute these here. Can potentially use log-sum-exp trick
        # as described here:
        # https://gregorygundersen.com/blog/2020/02/09/log-sum-exp/
        c_j = apply(X %*% beta[, -j, drop = FALSE], 1L, logsumexp)
        eta_j = X %*% beta[,j] - c_j
        for(i in 1:N){
          omega[i,j] <- BayesLogit::rpg(num=1, h = temper, z = eta_j[i]) # update omega
        }
        positive_finite(omega[, j], "Polya-Gamma draw")
        if (p<= N){ # Low-dimensional
          # Using Cholesky is more efficient than brute-force computation including direct matrix inverse using solve()
          # a term like "Sigma0_j %*% m0_j" is omitted because the prior mean m0_j term is zero.
          Z <- matrix(stats::rnorm(p), p, 1) # A p-vector of standard normals
          U <- chol(crossprod(x = X, y = (X * omega[,j])) + Sigma0_j) # Upper Cholesky factor of the inverse of V_j
          beta[,j] <- solve_chol(U, (t(X) %*% (delta[,j] + (omega[,j] * c_j)))) + backsolve(U, Z)
        }else if (p>N){# High-dimensional
          D = diag(lambda_sq[,j]*tau_sq, nrow = p, ncol = p)
          Phi = (X * sqrt(omega[,j]))
          alpha_s = sqrt(omega[,j]) * c_j + delta[,j] / sqrt(omega[,j])
          # Step 1.
          u = MASS::mvrnorm(n = 1, mu = rep(0, p), Sigma = D)
          e = MASS::mvrnorm(n = 1, mu = rep(0, N), Sigma = diag(N))
          # Step 2.
          v = Phi%*%u + e
          # Step 3.
          U = chol(crossprod(x = t(Phi), y = D%*%t(Phi)) + diag(N))
          w = solve_chol(U, (alpha_s-v)) # N by 1
          # Step 4.
          beta[,j] = u + D%*%t(Phi)%*%w
        }
      }
    }else if (prior_beta == "Gaussian"){
      # 1.1 Update beta and omega
      for(j in 1:(J-1)){
        c_j = apply(X %*% beta[, -j, drop = FALSE], 1L, logsumexp)
        eta_j = X %*% beta[,j] - c_j
        for(i in 1:N){
          omega[i,j] <- BayesLogit::rpg(num=1, h = temper, z = eta_j[i]) # update omega
        }
        positive_finite(omega[, j], "Polya-Gamma draw")
        if (p<=N){
          Z <- matrix(stats::rnorm(p), p, 1) # A p-vector of standard normals
          U <- chol(crossprod(x = X, y = (X * omega[,j])) + V0inv[,,j]) # Upper Cholesky factor of the inverse of V_j
          beta[,j] <- solve_chol(U, (t(X) %*% (delta[,j] + (omega[,j] * c_j)) + V0inv[,,j] %*% m0[,j])) + backsolve(U, Z)
        }else if (p>N){
          D = V0[,,j]
          Phi = (X * sqrt(omega[,j]))
          alpha_s = sqrt(omega[,j]) * c_j + delta[,j] / sqrt(omega[,j])
          # Step 1.
          u = MASS::mvrnorm(n = 1, mu = rep(0, p), Sigma = D)
          e = MASS::mvrnorm(n = 1, mu = rep(0, N), Sigma = diag(N))
          # Step 2.
          v = Phi%*%u + e
          # Step 3.
          U = chol(crossprod(x = t(Phi), y = D%*%t(Phi)) + diag(N))
          w = solve_chol(U, (alpha_s-v)) # N by 1
          # Step 4.
          beta[,j] = u + D%*%t(Phi)%*%w
        }
      }
    }
    # 3. Update kappa
    # Create the P matrix
    # P[i,j] is equal to p_{i,j} from the supplement of Polson et al. (2013).
    logP <- log_probabilities(X %*% beta)
    P <- exp(logP)
    if (is.character(prior_kappa) && prior_kappa == "default"){
      for(j in 1:(J-1)){
        # Now calculate the range of values that kappa[j] can take under its full
        # conditional distribution and the vector prob_vec, which is defined
        # such that prob_vec[k] is the log of the conditional probability that
        # kappa[j] = range[k] given all the other data and parameters.
        if(j==1){
          range <- candidate_range(0, kappa[j+1])
          lr <- length(range)
          log_prob_vec_prop <- rep(0, lr)
          for(l in 1:lr){
            log_prob_vec_prop[l] <- (
              temper*(sum(logP[1:range[l],j]) + sum(logP[(range[l]+1):kappa[j+1],j+1])) -
              length(1:range[l])*log(length(1:range[l])) -
              length((range[l]+1):kappa[j+1])*log(length((range[l]+1):kappa[j+1]))
            )
          }
          prob_vec <- normalize_weights(log_prob_vec_prop)
        } else if(j==L){
          range <- candidate_range(kappa[j-1], N)
          lr <- length(range)
          log_prob_vec_prop <- rep(0, lr)
          for(l in 1:lr){
            log_prob_vec_prop[l] <- (
              temper*(sum(logP[(kappa[j-1]+1):range[l],j]) + sum(logP[(range[l]+1):N,j+1])) -
              length((kappa[j-1]+1):range[l])*log(length((kappa[j-1]+1):range[l])) -
              length((range[l]+1):N)*log(length((range[l]+1):N))
            )
          }
          prob_vec <- normalize_weights(log_prob_vec_prop)
        } else{
          range <- candidate_range(kappa[j-1], kappa[j+1])
          lr <- length(range)
          log_prob_vec_prop <- rep(0, lr)
          for(l in 1:lr){
            log_prob_vec_prop[l] <- (
              temper*(sum(logP[(kappa[j-1]+1):range[l],j]) + sum(logP[(range[l]+1):kappa[j+1],j+1])) -
              length((kappa[j-1]+1):range[l])*log(length((kappa[j-1]+1):range[l])) -
              length((range[l]+1):kappa[j+1])*log(length((range[l]+1):kappa[j+1]))
            )
          }
          prob_vec <- normalize_weights(log_prob_vec_prop)
        }
        if(lr > 1){
          kappa[j] <- sample(x=range, size=1, prob=prob_vec) # update kappa
        } else {
          kappa[j] <- range[1L]
        }
      }
    }else if (is.function(prior_kappa)){ # else if the user have a joint PMF for the kappas
      for(j in 1:(J-1)){
        if(j==1){
          range <- candidate_range(0, kappa[j+1])
          lr <- length(range)
          log_prob_vec_prop <- rep(0, lr)
          for(l in 1:lr){
            kappa_prop = kappa
            kappa_prop[j] = range[l]
            jt_pmf_kappa = joint_prior_mass(kappa_prop)
            log_prob_vec_prop[l] <- temper*(sum(logP[1:range[l],j]) + sum(logP[(range[l]+1):kappa[j+1],j+1])) + log(jt_pmf_kappa)
          }
          prob_vec <- normalize_weights(log_prob_vec_prop)
        } else if(j==L){
          range <- candidate_range(kappa[j-1], N)
          lr <- length(range)
          log_prob_vec_prop <- rep(0, lr)
          for(l in 1:lr){
            kappa_prop = kappa
            kappa_prop[j] = range[l]
            jt_pmf_kappa = joint_prior_mass(kappa_prop)
            log_prob_vec_prop[l] <- temper*(sum(logP[(kappa[j-1]+1):range[l],j]) + sum(logP[(range[l]+1):N,j+1])) + log(jt_pmf_kappa)
          }
          prob_vec <- normalize_weights(log_prob_vec_prop)
        } else{
          range <- candidate_range(kappa[j-1], kappa[j+1])
          lr <- length(range)
          log_prob_vec_prop <- rep(0, lr)
          for(l in 1:lr){
            kappa_prop = kappa
            kappa_prop[j] = range[l]
            jt_pmf_kappa = joint_prior_mass(kappa_prop)
            log_prob_vec_prop[l] <- temper*(sum(logP[(kappa[j-1]+1):range[l],j]) + sum(logP[(range[l]+1):kappa[j+1],j+1])) + log(jt_pmf_kappa)
          }
          prob_vec <- normalize_weights(log_prob_vec_prop)
        }
        if(lr > 1){
          kappa[j] <- sample(x=range, size=1, prob=prob_vec) # update kappa
        } else {
          kappa[j] <- range[1L]
        }
      }
    } # End of the case if the user has a joint PMF of kappa
    # Update kappa in the scale of the entire series
    original_kappa = train_idx[kappa]
    #### Compute AUC for model selection ####
    if (model_selection  && (iter > num_warmup)){
      # "ratios" are the probabilities of a subject belonging to
      # class j+1 other than class j
      # Adjacent log odds have exactly the same ROC ranking as the
      # conditional probabilities, without saturation near zero or one.
      Y_all <- kappa_to_Y(original_kappa, n, L)
      ci = matrix(0, nrow = J-1, ncol = 3) # each row has a LB, point estimate, and UB
      for (j in 1:(J-1)){
        if(j == 1){
          rg = 1:original_kappa[j+1]
          lg = length(rg)
          ratios = as.numeric(X_all[rg, , drop = FALSE] %*% (beta[, j+1] - beta[, j]))
          Y_auc = Y_all[rg, j + 1]
          idx_test = which(rg %in% holdout_idx)
          P_auc = ratios[idx_test]
          Y_auc = Y_auc[idx_test]
          ci[j,] = auc_interval(Y_auc, P_auc)
        }else if (j == J-1){
          rg = (original_kappa[j-1]+1):n
          lg = length(rg)
          ratios = as.numeric(X_all[rg, , drop = FALSE] %*% (beta[, j+1] - beta[, j]))
          Y_auc = Y_all[rg, j + 1]
          idx_test = which(rg %in% holdout_idx)
          P_auc = ratios[idx_test]
          Y_auc = Y_auc[idx_test]
          ci[j,] = auc_interval(Y_auc, P_auc)
        }else{
          rg = (original_kappa[j-1]+1):original_kappa[j+1]
          lg = length(rg)
          ratios = as.numeric(X_all[rg, , drop = FALSE] %*% (beta[, j+1] - beta[, j]))
          Y_auc = Y_all[rg, j + 1]
          idx_test = which(rg %in% holdout_idx)
          P_auc = ratios[idx_test]
          Y_auc = Y_auc[idx_test]
          ci[j,] = auc_interval(Y_auc, P_auc)
        }
      }
    }
    if (prior_beta == "horseshoe"){
      if (temper == 1){
        if (model_selection){
          return( list(beta=beta, kappa=kappa, original_kappa=original_kappa, omega=omega, tau_sq = tau_sq, xi=xi, nu=nu, lambda_sq=lambda_sq, P=P, Auc=ci))
        }else{
          return( list(beta=beta, kappa=kappa, original_kappa=original_kappa, omega=omega, tau_sq = tau_sq, xi=xi, nu=nu, lambda_sq=lambda_sq, P=P) )
        }
      }else{
        return( list(beta=beta, kappa=kappa, original_kappa=original_kappa, omega=omega, tau_sq = tau_sq, xi=xi, nu=nu, lambda_sq=lambda_sq, P=P) )
      }
    }
    else if (prior_beta == "Gaussian"){
      if (model_selection){
        return( list(beta = beta, kappa = kappa, original_kappa = original_kappa, omega = omega, P = P, Auc = ci) )
      }else{
        return( list(beta = beta, kappa = kappa, original_kappa = original_kappa, omega = omega, P = P) )
      }
    }else{
      return( list(beta = beta, kappa = kappa, original_kappa = original_kappa, omega = omega) )
    }
  }

  #### Fitting single change point ####
  if (num_CP == 1){
    L = num_CP
    # priors for beta (logistic regression coefficients)
    if (prior_beta == "horseshoe"){
      xi = 1
      tau_sq = 1
      nu = rep(1, p)
      lambda_sq = rep(1, p)
    }else if (prior_beta == "Gaussian"){
      m0 <- rep(0, times = p) # m0[,j] is the prior mean vector
      sd0 <- sd_beta # prior variance
      V0 <- sd0*diag(p) # V0 is the prior cov. matrix
      V0inv <- 1/sd0*diag(p) # V0inv is the prior precision matrix
    }
    # Initialize working quantities
    beta = rep(0, times = p)
    omega = rep(0, times = N)
    P = rep(0, times = N)
    kappa = round(N/2) # Retain midpoint initialization for one changepoint.
    if (is.function(prior_kappa) && joint_prior_mass(kappa) <= 0)
      stop("The initial changepoint must have positive prior mass.")
    # Create an output list
    out = list()
    out$init_cp = kappa
    out$Beta = matrix(0, nrow = num_iter-num_warmup, ncol = p)
    out$P = matrix(0, nrow = num_iter-num_warmup, ncol = N)
    out$Kappa = rep(0, times = num_iter - num_warmup)
    out$AUC = rep(0, times = num_iter - num_warmup)
    # START of the for-loop #
    if (print_progress){
      pb = utils::txtProgressBar(min = 0, max = num_iter, style = 3)
      on.exit(close(pb), add = TRUE)
    }
    for (iter in 1:num_iter){
      # 1. Update omega
      eta = X %*% beta
      for (i in 1:N){
        omega[i] = if (tempering == 1) BayesLogit::rpg.devroye(num = 1, h = 1, z = eta[i]) else
          BayesLogit::rpg(num = 1, h = tempering, z = eta[i])
      }
      # Update delta
      positive_finite(omega, "Polya-Gamma draw")
      delta = tempering * c(rep(-0.5, times = kappa), rep(0.5, times = N-kappa))
      # 2. Update beta
      # Update HS prior parameters
      if (prior_beta =="horseshoe"){
        # Update scales in both low- and high-dimensional fits.
          xi <- LearnBayes::rigamma(n=1, a = 1, b = 1+1/tau_sq)
          xi_inv <- 1/xi
          tau_sq <- LearnBayes::rigamma(n=1, a = 0.5*(p+1), b =  xi_inv + 0.5*sum(beta^2 / lambda_sq))
          for (d in 1:p){
            nu[d] <- LearnBayes::rigamma(n=1, a = 1, b = 1+1/lambda_sq[d])
            lambda_sq[d] <- LearnBayes::rigamma(n=1, a = 1, b = 1/nu[d] + 0.5/tau_sq*beta[d]^2)
          }
        positive_finite(lambda_sq * tau_sq, "Horseshoe variance")
        if(p <= N){ # Low dimensional
          Sigma0 = diag(1/(lambda_sq*tau_sq), nrow = p, ncol = p)
          Z <- matrix(stats::rnorm(p), p, 1) # A p-vector of standard normals
          U <- chol(crossprod(x = X, y = (X * as.numeric(omega))) + Sigma0) # Upper Cholesky factor of the inverse of V_j
          beta <- solve_chol(U, (t(X)%*%delta)) + backsolve(U, Z)
        }else if (p>N){ # High dimensional
          Phi = (X * sqrt(as.numeric(omega)))
          D = diag(lambda_sq*tau_sq, nrow = p, ncol = p)
          alpha_s = delta / sqrt(as.numeric(omega))
          # Step 1.
          u = MASS::mvrnorm(n = 1, mu = rep(0, p), Sigma = D)
          e = MASS::mvrnorm(n = 1, mu = rep(0, N), Sigma = diag(N))
          # Step 2.
          v = Phi%*%u + e
          # Step 3.
          U = chol(crossprod(x = t(Phi), y = D%*%t(Phi)) + diag(N))
          w = solve_chol(U, (alpha_s-v)) # N by 1
          # Step 4.
          beta = u + D%*%t(Phi)%*%w
        }
      }else if (prior_beta == "Gaussian"){
        if (p <= N){ # Low dimensional
          Z <- matrix(stats::rnorm(p), p, 1) # A p-vector of standard normals
          U <- chol(crossprod(x = X, y = (X * as.numeric(omega))) + V0inv) # Upper Cholesky factor of the inverse of V_j
          beta <- solve_chol(U, (t(X)%*%delta + V0inv %*% m0)) + backsolve(U, Z)
        }else if(p>N){# High dimensional
          # Implement the fast sampler from https://arxiv.org/abs/1506.04778
          Phi = (X * sqrt(as.numeric(omega)))
          D = V0
          alpha_s = delta / sqrt(as.numeric(omega))
          # Step 1.
          u = MASS::mvrnorm(n = 1, mu = rep(0, p), Sigma = D)
          e = MASS::mvrnorm(n = 1, mu = rep(0, N), Sigma = diag(N))
          # Step 2.
          v = Phi%*%u + e
          # Step 3.
          U = chol(crossprod(x = t(Phi), y = D%*%t(Phi)) + diag(N))
          w = solve_chol(U, (alpha_s-v)) # N by 1
          # Step 4.
          beta = u + D%*%t(Phi)%*%w
        }
      }
      # 3.Update P and Kappa
      eta <- as.numeric(X %*% beta)
      if (any(!is.finite(eta))) stop("Non-finite binary linear predictor.")
      P <- stats::plogis(eta)
      logP1 <- stats::plogis(eta, log.p = TRUE)
      logP0 <- stats::plogis(eta, lower.tail = FALSE, log.p = TRUE)
      # range <- 1:(N-1)
      range <- candidate_range(0, N)
      lr <- length(range)
      log_prob_vec_prop <- rep(0, lr)
      if (is.character(prior_kappa) && prior_kappa == "default"){
        for(l in 1:lr){
          log_prob_vec_prop[l] = (
            tempering * (sum(logP0[1:range[l]]) + sum(logP1[(range[l]+1):N])) -
            length(1:range[l])*log(length(1:range[l])) -
            length((range[l]+1):N)*log(length((range[l]+1):N))
          )
        }
      }else if (is.function(prior_kappa)){
        jt_pmf_kappa = joint_prior_mass(kappa)
        for(l in 1:lr){
          jt_pmf_kappa = joint_prior_mass(range[l])
          log_prob_vec_prop[l] = tempering * (sum(logP0[1:range[l]]) + sum(logP1[(range[l]+1):N])) + log(jt_pmf_kappa)
        }
      }
      prob_vec <- normalize_weights(log_prob_vec_prop)
      if(lr > 1){
        kappa = sample(x=range, size=1, prob=prob_vec)
      } else {
        kappa = range[1L]
      }
      # 4. Compute AUC
      original_kappa = train_idx[kappa]
      # Update Y based on the kappa in the entire series
      Y = c(rep(0, times = original_kappa), rep(1, times = n - original_kappa))
      if (model_selection  && (iter > num_warmup)){
        Y_auc = Y[seq_len(length(Y)) %% thinning == 0]
        # Update the matrix of probabilities of success based on the fitted coefficients
        # Log odds retain the probability ranking without rounding to 0 or 1.
        P_auc <- as.numeric(X_auc %*% beta)
        auc_ci = auc_interval(Y_auc, P_auc)
      }
      # 5. Store samples after burn-in
      if(iter > num_warmup){
        out$Beta[iter-num_warmup,] = beta
        out$P[iter-num_warmup,] = P
        out$Kappa[iter-num_warmup] = original_kappa
        if (model_selection){
          out$AUC[iter-num_warmup] = auc_ci[1] # Take the CI lower bound
        }
      }
      if (print_progress){
        utils::setTxtProgressBar(pb, iter)
      }
    }
    end_time  = Sys.time()
    # Progress bar is closed on function exit, including on errors.
    # END of the for-loop #
    # Use posterior modes as changepoint estimates
    cp = table(out$Kappa)
    out$Kappa_mode = as.integer(names(cp)[which.max(cp)]) # posterior mode
    if (model_selection){
      out$init_cp = NULL
      out$P = NULL
      out$Beta = NULL
    }
    else{
      # Compute posterior mean coefficients
      out$Beta_mean = apply(out$Beta, 2, mean)
      out$P_mean = apply(out$P, 2, mean)
      out$AUC = NULL # no need for AUC if we are not doing model selection
      out$init_cp = NULL
    }
    if (print_outputs) utils::str(out)
    return(out)
  }

  #### Fitting multiple change points #####
  else if (num_CP >1){
    # Initialization before the big for-loop
    L = num_CP # the number of change points
    J = L+1 # the number of classes
    n = dim(data)[1] # the length of the entire series
    # priors for beta (multinomial logistic regression coefficients)
    if (prior_beta == "horseshoe"){
      xi = 1
      tau_sq = 1
      nu = matrix(1, nrow = p, ncol = J-1)
      lambda_sq = matrix(1, nrow = p, ncol = J-1)
    }else if (prior_beta == "Gaussian"){
      m0 <- array(data=0, dim=c(p,J-1)) # m0[,j] is the prior mean of the jth coef vec
      V0 <- array(data=0, dim=c(p,p,J-1)) # V0[,,j] is the prior cov. of the jth coef vec
      V0inv <- array(data=0, dim=c(p,p,J-1)) # V0inv[,,j] is the inverse of V0[,,j]
      sd0 <- sd_beta # prior variance
      for(j in 1:(J-1)){
        V0[,,j] <- sd0*diag(p)
        V0inv[,,j] <- 1/sd0*diag(p)
      }
    }
    beta= matrix(0, nrow = p, ncol = J)
    omega = matrix(0, nrow = N, ncol = J)
    P = matrix(0, nrow = N, ncol = J)
    if (init == "even"){
      kappa = rep(0, L)
      for (l in 1:L){
        kappa[l] = round(N/J) *l
      }
    }else if (init == "kcp"){
      algo = rpt$Dynp(model="rbf", min_size=as.integer(L))$fit(X)
      bkps =  algo$predict(n_bkps=as.integer(L))
      kappa = unlist(bkps[1:(length(bkps)-1)])
    }else if (init == "multirank"){
      algo = rpt$Dynp(model="rank", min_size=as.integer(L))$fit(X)
      bkps =  algo$predict(n_bkps=as.integer(L))
      kappa = unlist(bkps[1:(length(bkps)-1)])
    }else if (init == "ecp"){
      ecp_estimates = ecp::e.divisive(X, k = L)$estimates
      kappa = ecp_estimates[2:(L+1)] - 1L
    }
    if (!valid_initialization(kappa)) {
      if (init == "even") kappa <- floor(seq_len(L) * N / J)
      if (!valid_initialization(kappa))
        stop("Initialization violates the strict training-series support; use init = 'even'.")
    }
    if (is.function(prior_kappa) && joint_prior_mass(kappa) <= 0)
      stop("Initial changepoints must have positive prior mass.")
    out = list() # An output list
    out$init_cp = kappa # store the initial CP
    out$Kappa = matrix(data = rep(0, times = (num_iter-num_warmup)*L), nrow = num_iter-num_warmup, ncol = L)
    out$Beta = array(data = rep(x = 0, times = (num_iter-num_warmup)*p*J), dim = c((num_iter-num_warmup), p, J))
    out$P = array(data = rep(x = 0, times = (num_iter-num_warmup)*N*J), dim = c((num_iter-num_warmup), N, J))
    out$AUC = matrix(0, nrow = num_iter-num_warmup, ncol = L)
    # START of the for-loop fitting multiple change points #
    if (print_progress){
      pb = utils::txtProgressBar(min = 0, max = num_iter, style = 3)
      on.exit(close(pb), add = TRUE)
    }
    start = Sys.time()
    for (iter in 1:num_iter){
      # 1. Update samples
      gs = GS_update(temper = tempering, beta, kappa, omega, tau_sq, xi, nu, lambda_sq, X_all, X, train_idx, holdout_idx, n, N, p, L, J,V0inv, V0, model_selection)
      beta = gs$beta
      kappa = gs$kappa
      omega = gs$omega
      P = gs$P
      auc_ci = gs$Auc
      original_kappa = gs$original_kappa
      if (prior_beta == "horseshoe"){
        tau_sq = gs$tau_sq
        xi = gs$xi
        nu = gs$nu
        lambda_sq = gs$lambda_sq
      }
      # 2. Store the samples in the current iteration
      if ( iter > num_warmup){
        out$Beta[iter-num_warmup,,] = beta
        out$Kappa[iter-num_warmup,] = original_kappa
        out$P[iter-num_warmup,,] = P
        if (model_selection){
          out$AUC[iter-num_warmup,] = auc_ci[,1] # Take the CI lower bounds
        }
      }
      if (print_progress){ utils::setTxtProgressBar(pb, iter)}
    }
    end = Sys.time()
    # Progress bar is closed on function exit.
    # END of the for-loop fitting multiple change points #
    runTime = end - start
    # Use posterior modes as changepoint estimates
    cps = c()
    for (l in 1:L){
      cp = table(out$Kappa[,l])
      kappa = as.integer(names(cp)[which.max(cp)])
      cps = c(cps, kappa)
    }
    out$Kappa_mode = cps # posterior modes
    if (model_selection){
      out$init_cp = NULL
      out$P = NULL
      out$Beta = NULL
    }else{
      # Compute posterior mean coefficients
      out$Beta_mean = apply(out$Beta, c(2,3), mean)
      out$P_mean = apply(out$P, c(2,3), mean)
      out$AUC = NULL # no need for AUC if we are not doing model selection
      out$init_cp = NULL
    }
    if (print_outputs){
      cat("\n" ,"The algorithm took ", format(runTime), "to run ", num_iter,
          "iterations, including ", num_warmup, "burn-in iterations", "to fit ", num_CP,
          "change points in a data set in", p, "dimensions.", "\n" )
      utils::str(out)
    }
    return(out)
  }
}
