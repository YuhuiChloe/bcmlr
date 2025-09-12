bcmlr <- function(data, num_CP, init = "even", prior = "Gaussian", alpha_f = 0.1, 
                  min_size = 30, sd_beta = 3, tempering = 1, thinning = 1, 
                  num_iter = 5000, num_warmup = 2500, print_outputs = FALSE, 
                  print_progress = FALSE, model_selection = FALSE){
  #### Min size warning ####
  if (thinning >= min_size){
    stop("The minimum segment length (min_size) should be greater than the reciprocal of the held-out fraction (thinning).")
  }
  #### Convert min_size ####
  if (thinning != 1){
    min_size = round(min_size*(1 - 1/thinning))
  }
  #### NA/infinity warning ####
  if (anyNA(data) || any(is.infinite(unlist(data)))){
    stop("NAs or infinity exist in the input data.")
  }
  #### Check edge-cases ####
  if (length(thinning) > 2 || thinning < 1){
    stop("thinning must be a positive integer no less than 1")
  }
  if (init != "even" && init != "kcp" && init != "ecp" && init != "multirank"){
    stop("Initialization must be one of these: 'even', 'kcp', 'ecp', 'multirank'")
  }
  if (thinning > min_size){
    stop("Mininum segment length being smaller than thinning may cause ")
  }
  #### Center & standardize the data ####
  # A function to center & standardize a matrix x
  center_standardize <- function(x, m, sdv){
    for (i in 1:ncol(x)){
      x[,i] = x[,i] - m[i]
      x[,i] = x[,i]/sdv[i]
    }
    return(x)
  }
  X_all = as.matrix(data) 
  n = nrow(X_all) # the sample size of the input data
  #### Hold out samples ####
  if (thinning != 1){
    all_idx = seq_len(n) # all indices from the input data
    holdout_idx = all_idx[all_idx%%thinning == 0] # the indices of the held-out samples
    train_idx = all_idx[-holdout_idx] # the indices for fitting the model (the trainng set indices)
    X_auc = X_all[holdout_idx, , drop = FALSE] # rows of X selected to calculate AUC
    X = X_all[-holdout_idx, , drop = FALSE] # rows of X selected to fit the model  
    if (ncol(X) == 1){
      m = mean(X, na.rm = TRUE)
    }else{
      m = colMeans(X, na.rm = TRUE) # column means in the training set
    }
    sdv = apply(X, 2, FUN = sd) # standard deviation of the training set
    sdv[sdv == 0] = 1 # Avoid divide-by-zero in standardization: sdv = 0 could occur if the input data has a column of the same constant. 
    X = center_standardize(X, m, sdv) # center & standardize the subseries used to fit the model with its column means
    X_auc = center_standardize(X_auc, m, sdv) # center & standardize the held-out samples with the same column means
    # Update the original data set with the centered-standardized X and X_auc
    X_all[holdout_idx,] = X_auc
    X_all[train_idx,] = X
    # Quick sanity checks: make sure all the indices match 
    stopifnot(
      nrow(X_auc) == length(holdout_idx),
      nrow(X)     == length(train_idx),
      isTRUE(all.equal(X_all[holdout_idx, ], X_auc)),
      isTRUE(all.equal(X_all[train_idx, ], X))
    )
  }else{
    all_idx = seq_len(n)
    holdout_idx = all_idx 
    train_idx = all_idx
    if (dim(X_all)[2] == 1){
      m = mean(X_all, na.rm = TRUE)
    }else{
      m = colMeans(X_all, na.rm = TRUE) # column means in the training set
    }
    sdv = apply(X_all, 2, FUN = sd)
    sdv[sdv == 0] = 1 # Avoid divide-by-zero in standardization: sdv = 0 could occur if the input data has a column of the same constant. 
    X = center_standardize(X_all, m, sdv)
    X_all = X
    # No need for X_auc since we are not computing AUCs for model selection
  }
  N = dim(X)[1] # Sample size of training set
  p = dim(X)[2] # Dimension of the training set
  #### Utility Functions ####
  # https://arxiv.org/pdf/2205.04997
  # https://gregorygundersen.com/blog/2020/02/09/log-sum-exp/
  # A function that mitigates computational issues in massive division or multiplication of extreme values.
  logsumexp <- function(x){
    c <- max(x)
    return(c + log(sum(exp(x - c))))
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
    Y <- model.matrix(~ -1 + class, data=class)
    return(Y)
  }
  # A function that updates \beta, \kappa, and  \omega per iteration  #
  GS_update <- function(temper, beta, kappa, omega, tau_sq, xi, nu, lambda, X_all, X, train_idx, holdout_idx, n, N, p, L, J, model_selection){
    # 0. Update Y, delta 
    Y = kappa_to_Y(kappa, N, L)
    delta = temper*(Y-1/2)
    # 1. Update HS prior (on beta) parameters & beta & omega
    if (prior == "horseshoe"){
      # 1.1 Update HS priors
      xi_inv <- rgamma(n=1, shape = 1, scale = 1+1/tau_sq)
      xi <- 1/xi_inv
      for (j in 1:(J-1)){
        tau_sq <- 1/rgamma(n=1, shape = 0.5*(p+1), scale =  xi_inv + 0.5*sum(beta[,j]^2 / lambda[,j]^2))
        for (d in 1:p){
          nu[d,j] <- 1/rgamma(n=1, shape = 1, scale = 1+1/lambda[d,j]^2)
          lambda[d,j] <- 1/rgamma(n=1, shape = 1, scale = 1/nu[d,j] + 0.5/tau_sq*beta[d,j]^2)
        }
      }
      # 1.2 Update beta and omega 
      for(j in 1:(J-1)){
        Sigma0_j = diag(lambda[,j]^2*tau_sq, nrow = p, ncol = p) 
        # Polson et al. (2013) describes or implies the existence of matrices
        # C and eta. We compute these here. Can potentially use log-sum-exp trick
        # as described here:
        # https://gregorygundersen.com/blog/2020/02/09/log-sum-exp/
        c_j = log(rowSums(exp(X %*% beta[,-j])))
        eta_j = X %*% beta[,j] - c_j
        for(i in 1:N){
          omega[i,j] <- rpg(num=1, h = temper, z = eta_j[i]) # update omega
        }
        # Using Cholesky is more efficient than brute-force computation including direct matrix inverse using solve()
        # a term like "Sigma0_j %*% m0_j" is omitted because the prior mean term is zero. 
        Z <- matrix(rnorm(p), p, 1) # A p-vector of standard normals
        U <- chol(crossprod(x = X, y = diag(omega[,j]) %*% X) + Sigma0_j) # Upper Cholesky factor of the inverse of V_j
        beta[,j] <- chol2inv(U) %*% (t(X) %*% (delta[,j] + diag(omega[,j]) %*% c_j)) + backsolve(U, Z)
      }
    }else if (prior == "Gaussian"){
      # 1.1 Update beta and omega 
      for(j in 1:(J-1)){
        c_j = log(rowSums(exp(X %*% beta[,-j])))
        eta_j = X %*% beta[,j] - c_j
        for(i in 1:N){
          omega[i,j] <- rpg(num=1, h = temper, z = eta_j[i]) # update omega
        }
        Z <- matrix(rnorm(p), p, 1) # A p-vector of standard normals
        U <- chol(crossprod(x = X, y = diag(omega[,j]) %*% X) + V0inv[,,j]) # Upper Cholesky factor of the inverse of V_j
        beta[,j] <- chol2inv(U) %*% (t(X) %*% (delta[,j] + diag(omega[,j]) %*% c_j) + V0inv[,,j] %*% m0[,j]) + backsolve(U, Z)
      }
    }
    # 3. Update kappa (and Y, which is a function of kappa) 
    # Create the P matrix
    # P[i,j] is equal to p_{i,j} from the supplement of Polson et al. (2013).
    Phi_exp <- exp(X %*% beta) # Element-wise exponent of Phi matrix from supplement
    P <- t(apply(Phi_exp, 1, function(x) x/sum(x))) # Normalize the rows of Phi_exp to get P
    for(j in 1:(J-1)){
      # Now calculate the range of values that kappa[j] can take under its full
      # conditional distribution and the vector prob_vec, which is defined
      # such that prob_vec[k] is the log of the conditional probability that
      # kappa[j] = range[k] given all the other data and parameters.
      if(j==1){
        #range <- 1:(kappa[j+1]-1)
        range <- (1 + min_size):(kappa[j+1] - 1 - min_size)
        lr <- length(range)
        log_prob_vec_prop <- rep(0, lr)
        for(l in 1:lr){
          log_prob_vec_prop[l] <- temper*(sum(log(P[1:range[l],j])) + sum(log(P[(range[l]+1):kappa[j+1],j+1])))
          - length(1:range[l])*log(length(1:range[l]))
          - length((range[l]+1):kappa[j+1])*log(length((range[l]+1):kappa[j+1]))
        }
        #prob_vec <- exp(log_prob_vec_prop)/sum(exp(log_prob_vec_prop))
        prob_vec <- exp(log_prob_vec_prop - logsumexp(log_prob_vec_prop))
      } else if(j==L){
        #range <- (kappa[j-1]+1):(N-1)
        range <- (kappa[j-1]+1 + min_size):(N - 1 - min_size)
        lr <- length(range)
        log_prob_vec_prop <- rep(0, lr)
        for(l in 1:lr){
          log_prob_vec_prop[l] <- temper*(sum(log(P[(kappa[j-1]+1):range[l],j])) + sum(log(P[(range[l]+1):N,j+1])))
          - length((kappa[j-1]+1):range[l])*log(length((kappa[j-1]+1):range[l]))
          - length((range[l]+1):N)*log(length((range[l]+1):N))
        }
        #prob_vec <- exp(log_prob_vec_prop)/sum(exp(log_prob_vec_prop))
        prob_vec <- exp(log_prob_vec_prop - logsumexp(log_prob_vec_prop))
      } else{
        #range <- (kappa[j-1]+1):(kappa[j+1]-1)
        range <- (kappa[j-1]+1 + min_size):(kappa[j+1]-1 - min_size)
        lr <- length(range)
        log_prob_vec_prop <- rep(0, lr)
        for(l in 1:lr){
          log_prob_vec_prop[l] <- temper*(sum(log(P[(kappa[j-1]+1):range[l],j])) + sum(log(P[(range[l]+1):kappa[j+1],j+1])))
          - length((kappa[j-1]+1):range[l])*log(length((kappa[j-1]+1):range[l]))
          - length((range[l]+1):kappa[j+1])*log(length((range[l]+1):kappa[j+1]))
        }
        #prob_vec <- exp(log_prob_vec_prop)/sum(exp(log_prob_vec_prop))
        prob_vec <- exp(log_prob_vec_prop - logsumexp(log_prob_vec_prop))
      }
      if(lr > 1){
        kappa[j] <- sample(x=range, size=1, prob=prob_vec) # update kappa
      }
    }
    # Update kappa in the scale of the entire series
    original_kappa = train_idx[kappa] 
    #### Compute AUC for model selection ####
    if (model_selection){
      # "ratios" are the probabilities of a subject belonging to
      # class j+1 other than class j
      Phi_exp_all <- exp(X_all %*% beta) # Element-wise exponent of Phi matrix from supplement
      P_all <- t(apply(Phi_exp_all, 1, function(x) x/sum(x))) # Normalize the rows of Phi_exp to get P
      Y_all <- kappa_to_Y(original_kappa, n, L)
      ci = matrix(0, nrow = J-1, ncol = 3) # each row has a LB, point estimate, and UB
      for (j in 1:(J-1)){
        if(j == 1){
          rg = 1:original_kappa[j+1]
          lg = length(rg)
          ratios = P_all[rg,j+1]/(P_all[rg,j] + P_all[rg,j+1])
          Y_auc = Y_all[rg, j + 1]
          idx_test = which(rg %in% holdout_idx)
          P_auc = ratios[idx_test]
          Y_auc = Y_auc[idx_test]
          ci[j,] = as.numeric(suppressWarnings(ci.auc(as.factor(Y_auc), as.numeric(P_auc), conf.level= 1-alpha_f, quiet = TRUE)))
        }else if (j == J-1){
          rg = (original_kappa[j-1]+1):n
          lg = length(rg)
          ratios = P_all[rg,j+1]/(P_all[rg,j] + P_all[rg, j+1])
          Y_auc = Y_all[rg, j + 1]
          idx_test = which(rg %in% holdout_idx)
          P_auc = ratios[idx_test]
          Y_auc = Y_auc[idx_test]
          ci[j,] = as.numeric(suppressWarnings(ci.auc(as.factor(Y_auc), as.numeric(P_auc), conf.level= 1-alpha_f, quiet = TRUE)))
        }else{
          rg = (original_kappa[j-1]+1):original_kappa[j+1]
          lg = length(rg)
          ratios = P_all[rg,j+1]/(P_all[rg,j] + P_all[rg, j+1])
          Y_auc = Y_all[rg, j + 1]
          idx_test = which(rg %in% holdout_idx)
          P_auc = ratios[idx_test]
          Y_auc = Y_auc[idx_test]
          ci[j,] = as.numeric(suppressWarnings(ci.auc(as.factor(Y_auc), as.numeric(P_auc), conf.level= 1-alpha_f, quiet = TRUE)))
        }
      }
    }
    
    if (prior == "horseshoe"){
      if (temper == 1){
        if (model_selection){
          return( list(beta=beta, kappa=kappa, original_kappa=original_kappa, omega=omega, tau_sq = tau_sq, xi=xi, nu=nu, lambda=lambda, P=P, Auc=ci))
        }else{
          return( list(beta=beta, kappa=kappa, original_kappa=original_kappa, omega=omega, tau_sq = tau_sq, xi=xi, nu=nu, lambda=lambda, P=P) )
        }
      }else{
        return( list(beta=beta, kappa=kappa, original_kappa=original_kappa, omega=omega, tau_sq = tau_sq, xi=xi, nu=nu, lambda=lambda) )
      }
    }
    else if (prior == "Gaussian"){
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
    if (prior == "horseshoe"){
      xi = 1
      tau_sq = 1
      nu = rep(1, p)
      lambda = rep(1, p)
    }else if (prior == "Gaussian"){
      m0 <- rep(0, times = p) # m0[,j] is the prior mean vector
      sd0 <- sd_beta # prior standard deviation
      V0 <- sd0*diag(p) # V0 is the prior cov. matrix
      V0inv <- 1/sd0*diag(p) # V0inv is the prior precision matrix
    }
    # Initialize working quantities
    beta = rep(0, times = p) 
    omega = rep(0, times = N) 
    P = rep(0, times = N) 
    kappa = round(N/2) # bclr does not have the multi-modality challenge. Use even initialization.  
    # Create an output list
    out = list()
    out$init_cp = kappa
    out$Beta = matrix(0, nrow = num_iter-num_warmup, ncol = p)
    out$P = matrix(0, nrow = num_iter-num_warmup, ncol = N)
    out$Kappa = rep(0, times = num_iter - num_warmup)
    out$AUC = rep(0, times = num_iter - num_warmup)
    # START of the for-loop #
    start_time = Sys.time()
    if (print_progress){
      pb = txtProgressBar(min = 0, max = num_iter, style = 3)
    }
    for (iter in 1:num_iter){
      # 1. Update omega 
      eta = X %*% beta
      for (i in 1:N){
        omega[i] = rpg.devroye(num = 1, h = 1, z = eta[i]) 
      }
      # 2. Update beta 
      # Update HS prior parameters 
      if (prior =="horseshoe"){
        xi_inv <- rgamma(n=1, shape = 1, scale = 1+1/tau_sq)
        xi <- 1/xi_inv
        tau_sq <- 1/rgamma(n=1, shape = 0.5*(p+1), scale =  xi_inv + 0.5*sum(beta^2 / lambda^2))
        for (d in 1:p){
          nu[d] <- 1/rgamma(n=1, shape = 1, scale = 1+1/lambda[d]^2)
          lambda[d] <- 1/rgamma(n=1, shape = 1, scale = 1/nu[d] + 0.5/tau_sq*beta[d]^2)
        }
        Sigma0 = diag(lambda^2*tau_sq, nrow = p, ncol = p)
        delta = c(rep(-0.5, times = kappa), rep(0.5, times = N-kappa))
        Z <- matrix(rnorm(p), p, 1) # A p-vector of standard normals
        U <- chol(crossprod(x = X, y = diag(omega) %*% X) + Sigma0) # Upper Cholesky factor of the inverse of V_j
        beta <- chol2inv(U) %*% (t(X)%*%delta) + backsolve(U, Z)
      }else if (prior == "Gaussian"){
        delta = c(rep(-0.5, times = kappa), rep(0.5, times = N-kappa)) 
        Z <- matrix(rnorm(p), p, 1) # A p-vector of standard normals
        U <- chol(crossprod(x = X, y = diag(omega) %*% X) + V0inv) # Upper Cholesky factor of the inverse of V_j 
        beta <- chol2inv(U) %*% (t(X)%*%delta + V0inv %*% m0) + backsolve(U, Z)
      }
      # 3.Update P and Kappa
      Phi_exp <- exp(X %*% beta) 
      P = Phi_exp/(1+Phi_exp)
      # range <- 1:(N-1)
      range <- (1+min_size):(N-1-min_size)
      lr <- length(range)
      log_prob_vec_prop <- rep(0, lr) 
      for(l in 1:lr){
        log_prob_vec_prop[l] = sum(log(1-P[1:range[l]])) + sum(log(P[(range[l]+1):N]))
        - length(1:range[l])*log(length(1:range[l]))
        - length((range[l]+1):N)*log(length((range[l]+1):N))
      } 
      prob_vec <- exp(log_prob_vec_prop - logsumexp(log_prob_vec_prop)) 
      if(lr > 1){
        kappa = sample(x=range, size=1, prob=prob_vec)
      }
      # 4. Compute AUC
      original_kappa = train_idx[kappa]
      # Update Y based on the kappa in the entire series
      Y = c(rep(0, times = original_kappa), rep(1, times = n - original_kappa))
      if (model_selection){
        Y_auc = Y[seq_len(length(Y)) %% thinning == 0]
        # Update the matrix of probabilities of success based on the fitted coefficients
        Phi_exp_auc <- exp(X_auc %*% beta) 
        P_auc = Phi_exp_auc/(1+Phi_exp_auc)
        #auc = suppressWarnings(ci.auc(as.factor(Y_auc), as.numeric(P_auc), conf.level= 1-alpha_f, quiet = TRUE)[1])
        auc_ci = as.numeric(suppressWarnings(ci.auc(as.factor(Y_auc), as.numeric(P_auc), conf.level= 1-alpha_f, quiet = TRUE)))
      }
      # 5. Store samples after burn-in
      if(iter > num_warmup){
        out$Beta[iter-num_warmup,] = beta 
        out$P[iter-num_warmup,] = P
        out$Kappa[iter-num_warmup] = original_kappa
        if (model_selection){
          out$AUC[iter-num_warmup] = auc_ci[1]
        }
      }
      if (print_progress){
        setTxtProgressBar(pb, iter)
      }
    }
    end_time  = Sys.time()
    if (print_progress){
      close(pb) 
    }
    # END of the for-loop #
    # Use posterior modes as changepoint estimates
    cp = table(out$Kappa)
    out$kappa_mode = as.integer(names(cp)[which.max(cp)]) # posterior mode 
    if (model_selection){
      out$init_cp = NULL
      out$P = NULL
      out$Beta = NULL
    }
    else{
      # Compute posterior mean coefficients
      out$Beta = apply(out$Beta, 2, mean)
      out$P = apply(out$P, 2, mean)
      out$AUC = NULL # no need for AUC if we are not doing model selection
      out$init_cp = NULL
    }
    return(out)
  }
  
  #### Fitting multiple change points #####
  else if (num_CP >1){
    # Initialization before the big for-loop 
    L = num_CP # the number of change points
    J = L+1 # the number of classes
    n = dim(data)[1] # the length of the entire series
    # priors for beta (multinomial logistic regression coefficients)
    if (prior == "horseshoe"){
      xi = 1
      tau_sq = 1
      nu = matrix(1, nrow = p, ncol = J-1)
      lambda = matrix(1, nrow = p, ncol = J-1)
    }else if (prior == "Gaussian"){
      m0 <- array(data=0, dim=c(p,J-1)) # m0[,j] is the prior mean of the jth coef vec
      V0 <- array(data=0, dim=c(p,p,J-1)) # V0[,,j] is the prior cov. of the jth coef vec
      V0inv <- array(data=0, dim=c(p,p,J-1)) # V0inv[,,j] is the inverse of V0[,,j]
      sd0 <- sd_beta # prior standard deviation
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
      algo = rpt$Dynp(model="rbf", min_size=L)$fit(X)
      bkps =  algo$predict(n_bkps=L)  
      kappa = unlist(bkps[1:(length(bkps)-1)])
    }
    else if (init == "multirank"){
      algo = rpt$Dynp(model="rank", min_size=L)$fit(X)
      bkps =  algo$predict(n_bkps=L)
      kappa = unlist(bkps[1:(length(bkps)-1)])
    }else if (init == "ecp"){
      ecp_estimates = e.divisive(X, k = L)$estimates
      kappa = ecp_estimates[2:(L+1)]
    }
    out = list() # An output list
    out$init_cp = kappa # store the initial CP
    out$Kappa = matrix(data = rep(0, times = (num_iter-num_warmup)*L), nrow = num_iter-num_warmup, ncol = L)
    out$Beta = array(data = rep(x = 0, times = (num_iter-num_warmup)*p*J), dim = c((num_iter-num_warmup), p, J))
    out$P = array(data = rep(x = 0, times = (num_iter-num_warmup)*N*J), dim = c((num_iter-num_warmup), N, J))
    out$AUC = matrix(0, nrow = num_iter-num_warmup, ncol = L)
    # START of the for-loop fitting multiple change points #
    if (print_progress){
      pb = txtProgressBar(min = 0, max = num_iter, style = 3)
    }
    start = Sys.time()
    for (iter in 1:num_iter){
      # 1. Update samples
      gs = GS_update(temper = tempering, beta, kappa, omega, tau_sq, xi, nu, lambda, X_all, X, train_idx, holdout_idx, n, N, p, L, J, model_selection)
      beta = gs$beta
      kappa = gs$kappa
      omega = gs$omega
      P = gs$P
      auc_ci = gs$Auc
      original_kappa = gs$original_kappa
      if (prior == "horseshoe"){
        tau_sq = gs$tau_sq
        xi = gs$xi
        nu = gs$nu 
        lambda = gs$lambda
      }
      # 2. Store the samples in the current iteration
      if ( iter > num_warmup){
        out$Beta[iter-num_warmup,,] = beta
        out$Kappa[iter-num_warmup,] = original_kappa
        out$P[iter-num_warmup,,] = P
        if (model_selection){
          out$AUC[iter-num_warmup,] = auc_ci[,1] 
        }
      }
      if (print_progress){ setTxtProgressBar(pb, iter)}
    }
    end = Sys.time()
    if (print_progress){ close(pb) }
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
      out$Beta = apply(out$Beta, c(2,3), mean)
      out$P = apply(out$P, c(2,3), mean)
      out$AUC = NULL # no need for AUC if we are not doing model selection
      out$init_cp = NULL
    }
    if (print_outputs){
      cat("\n" ,"The algorithm took ", format(runTime), "to run ", num_iter,
          "iterations, including ", num_warmup, "burn-in iterations", "to fit ", num_CP,
          "change points in a data set in", p, "dimensions.", "\n" )
      cat(str(out))
    }
    return(out)
  }
}
