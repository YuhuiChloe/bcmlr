bcmlr <- function(data, num_iter = 10000, num_warmup = 5000, init_CP = NA, num_CP,thinning, print_outputs = TRUE){

  X = as.matrix(data)
  if (anyNA(init_CP) == FALSE && length(init_CP) != num_CP){
    print("initial changepoints must have the same length with the number of changepoints")
  }
  
  if (length(thinning) >2 || thinning < 1){
    print("thinning must be a positive integer no less than 1")
    break
  }
  
  if (thinning != 1){
    X_auc = X[seq_len(nrow(X)) %% thinning == 0, ] # rows of X selected to calculate AUC
    X_auc = scale(X_auc, center = TRUE, scale = FALSE)
    X = X[seq_len(nrow(X)) %% thinning != 0, ] # rows of X selected to fit bclr 
    X = scale(X, center = TRUE, scale = FALSE)
  }else{
    X_auc = scale(X, center = TRUE, scale = FALSE)  # rows of X selected to calculate AUC
    X = scale(X, center = TRUE, scale = FALSE) # rows of X selected to fit bclr 
  }
  
  N = dim(X)[1] # Sample size
  p = dim(X)[2] # data dimension
  
  
  ### Functions used for both cases (single CP and multiple CPs) ####
  
  # https://arxiv.org/pdf/2205.04997
  # https://gregorygundersen.com/blog/2020/02/09/log-sum-exp/
  # A function that mitigate computational issues in 
  # massive division or multiplication of extreme values. 
  logsumexp <- function(x){
    c <- max(x)
    return(c + log(sum(exp(x - c))))
  }
  
  # A function that inputs binary labels and predicted probabilities
  # and returns an AUC value 
  calculate_roc_auc <- function(Y, P) {
    
    # Y: labels, P: predicted probabilities
    df <- data.frame(Y = Y, P = P)
    df <- df[order(-df$P),]
    
    # Compute the total positives and negatives
    total_positives <- sum(df$Y == 1)
    total_negatives <- sum(df$Y == 0)
    # Initialize variables to store the true positive rate (TPR) and false positive rate (FPR)
    TPR <- rep(0, times = total_positives)
    FPR <- rep(0, times = total_negatives)
    TP <- 0
    FP <- 0
    
    # Iterate over each threshold (sorted predicted probabilities)
    for (i in 1:nrow(df)) {
      if (df$Y[i] == 1) {
        TP <- TP + 1
      } else {
        FP <- FP + 1
      }
      # Calculate TPR and FPR for each threshold
      TPR[i] <- TP / total_positives
      FPR[i] <- FP / total_negatives
    }
    # Add the (0,0) and (1,1) points to the ROC curve
    TPR <- c(0, TPR, 1)
    FPR <- c(0, FPR, 1)
    
    # Calculate AUC using trapezoidal rule
    auc <- sum(diff(FPR) * (head(TPR, -1) + tail(TPR, -1)) / 2)
    
    return(auc)
  }

  #### Utility functions and functions for Gibbs sampling, Metropolis-Hastings check, and parallel computing ####
  
  # 1-1 Correspondence
  # A function that converts a change point vector to a response matrix with 
  # the change point structure
  kappa_to_Y <- function(kappa, N, L){
    # L <- length(kappa)
    class <- rep(1, kappa[1])
    for(j in 2:L){
      class <- c(class, rep(j, kappa[j] - kappa[j-1]))
    }
    class <- c(class, rep(L+1, N - kappa[L]))
    class <- as.factor(class)
    Y <- model.matrix(~ -1 + class, data=class)
    return(Y)
  } 
  
  # Use this function when calculating the log target in MH steps. 
  # A function that converts kappa (with length L) to L+1 penalization ratios. 
  kappa_to_alpha <- function(kappa, N, L){
    alpha = rep(0, times = L+1)
    for (j in 1:(L+1)){
      if (j == 1){
        alpha[j] = N/kappa[1]
      }else if(j == L+1){
        alpha[j] = N/(N-kappa[L])
      }else{
        alpha[j] = N/(kappa[j] - kappa[j-1])
      }
    }
    return(alpha)
  }
  
  
  #### A function that updates \beta, \kappa, and  \omega per iteration  ########
  GS_update <- function(temper, beta, kappa, omega, data, X, X_auc, N, L, J, m0, V0inv){
    
    # update Y, delta
    Y = kappa_to_Y(kappa, N, L)
    delta = temper*(Y-1/2)
  
    ### Update beta and omega #####
    for(j in 1:(J-1)){
      
      # Polson et al. (2013) describes or implies the existence of matrices 
      # C and eta. We compute these here. Can potentially use log-sum-exp trick
      # as described here: 
      # https://gregorygundersen.com/blog/2020/02/09/log-sum-exp/
      c_j = log(rowSums(exp(X %*% beta[,-j])))
      eta_j = X %*% beta[,j] - c_j 
      
      for(i in 1:N){ 
        omega[i,j] <- rpg(num=1, h = temper, z = eta_j[i]) # update omega
      }
      
      # V_j <- solve(t(X) %*% diag(omega[,j]) %*% X + V0inv[,,j])
      # V_j <- chol2inv(chol(t(X) %*% diag(omega[,j]) %*% X + V0inv[,,j]))
      # m_j <- V_j %*% (t(X) %*% (delta[,j] + diag(omega[,j]) %*% c_j) + V0inv[,,j] %*% m0[,j])
      # beta[,j] = mvrnorm(n = 1, mu = m_j, Sigma = V_j)  # update beta
      Z <- matrix(rnorm(p), p, 1) # A p-vector of standard normals
      U <- chol(crossprod(x = X, y = diag(omega[,j]) %*% X) + V0inv[,,j]) # Upper Cholesky factor of the inverse of V_j 
      beta[,j] <- chol2inv(U) %*% (t(X) %*% (delta[,j] + diag(omega[,j]) %*% c_j) + V0inv[,,j] %*% m0[,j]) + backsolve(U, Z)
    }

    
    #### Update kappa (and Y, which is a function of kappa) #### 
    
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
        range <- 1:(kappa[j+1]-1)
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
        range <- (kappa[j-1]+1):(N-1)
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
        range <- (kappa[j-1]+1):(kappa[j+1]-1)
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
    
    # "ratios" are the probabilities of a subject belonging to 
    # class j+1 other than class j
    X_all = as.matrix(data)
    Phi_exp_all <- exp(X_all %*% beta) # Element-wise exponent of Phi matrix from supplement
    P_all <- t(apply(Phi_exp_all, 1, function(x) x/sum(x))) # Normalize the rows of Phi_exp to get P
    
    n = dim(data)[1]
    idx_all = 1:n
    if (thinning ==1){
      idx_fit = idx_all
    }else{
      idx_fit = idx_all[idx_all%%thinning != 0]
    }
    original_kappa = idx_fit[kappa]
    Y_all <- kappa_to_Y(original_kappa, n, L)
    auc_values = rep(0, times = J-1)
    
    for (j in 1:(J-1)){
      if(j == 1){
        rg = 1:(original_kappa[j+1]-1)
        lg = length(rg)
        ratios = P_all[rg,j+1]/(P_all[rg,j] + P_all[rg,j+1])
        ratios = ratios[seq_len(lg) %% thinning == 0]
        Y_auc = Y_all[rg, j + 1]
        Y_auc = Y_auc[seq_len(lg) %% thinning == 0]
        auc_values[j] = calculate_roc_auc(Y_auc, ratios)
      }else if (j == J-1){
        rg = (original_kappa[j-1]+1):(n-1)
        lg = length(rg)
        ratios = P_all[rg,j+1]/(P_all[rg,j] + P_all[rg, j+1])
        ratios = ratios[seq_len(lg) %% thinning == 0]
        Y_auc = Y_all[rg, j + 1]
        Y_auc = Y_auc[seq_len(lg) %% thinning == 0]
        auc_values[j] = calculate_roc_auc(Y_auc, ratios)
      }else{
        rg = (original_kappa[j-1]+1):(original_kappa[j+1]-1)
        lg = length(rg)
        ratios = P_all[rg,j+1]/(P_all[rg,j] + P_all[rg, j+1])
        ratios = ratios[seq_len(lg) %% thinning == 0]
        Y_auc = Y_all[rg, j + 1]
        Y_auc = Y_auc[seq_len(lg) %% thinning == 0]
        auc_values[j] = calculate_roc_auc(Y_auc, ratios)
      }
    }
    
    if (temper == 1){
      return( list(beta = beta, kappa = kappa, original_kappa = original_kappa, omega = omega, P = P, Auc = auc_values) )
    }else{
      return( list(beta = beta, kappa = kappa, original_kappa = original_kappa, omega = omega) )
    }
    
  }
  
  ##### A function that return a value of the logarithm of the target density kernel ######
  logtarget <- function(temper, beta, kappa, omega, X, N, L, J, m0, V0inv){
    
    # update Y, delta, and alpha
    Y <- kappa_to_Y(kappa, N, L)
    delta <- temper*(Y-1/2)
    alpha <- kappa_to_alpha(kappa, N,L)
    
    logkernel = 0 
    for (j in 1:(J-1)){
      c_j = log(rowSums(exp(X %*% beta[,-j])))
      eta_j = X %*% beta[,j] - c_j 
      z_j = delta[,j]/omega[,j]
      
      prior_j = - 0.5*t(beta[,j]-m0[,j])%*%V0inv[,,j]%*%(beta[,j]-m0[,j]) 
      like_j = - 0.5*t(X%*%beta[,j] -c_j - z_j)%*%diag(omega[,j])%*%(X%*%beta[,j] -c_j - z_j) + 0.5*t(z_j)%*%delta[,j]
      logkernel = logkernel + prior_j + like_j
      
      for (i in 1:N){ 
        logkernel = logkernel + polyagamma_logpdf(x = omega[i,j], h = temper, z = 0) 
      }
    }
    
    penalty = temper*kappa[1]*log(alpha[1]) + temper*(N-kappa[L])*log(alpha[L+1])
    for(l in 2:L){
      penalty = penalty + temper*(kappa[l]-kappa[l-1])*log(alpha[l])
    }
    logkernel =  logkernel + penalty
    
    return(logkernel)
    
  }
  
  
  ##################################### Fitting multiple change points ####################################
  cat("\n", "Fitting ", num_CP ," change points......")
  
  L = num_CP # the number of change points 
  J = L+1 # the number of classes
  n = dim(data)[1] # the length of the entire series
  
  # priors for beta (multinomial logistic regression coefficients)
  m0 <- array(data=0, dim=c(p,J-1)) # m0[,j] is the prior mean of the jth coef vec
  V0 <- array(data=0, dim=c(p,p,J-1)) # V0[,,j] is the prior cov. of the jth coef vec
  V0inv <- array(data=0, dim=c(p,p,J-1)) # V0inv[,,j] is the inverse of V0[,,j]
  sd0 <- 3 # prior standard deviation 
  for(j in 1:(J-1)){
    V0[,,j] <- sd0*diag(p)
    V0inv[,,j] <- 1/sd0*diag(p)
  }
  
  
  # Initialize working quantities
  beta= matrix(0, nrow = p, ncol = J)
  omega = matrix(0, nrow = N, ncol = J)
  P = matrix(0, nrow = N, ncol = J)
  if (anyNA(init_CP)){
    kappa = rep(0, L)
    for (l in 1:L){
      kappa[l] = round(N/J) *l
    }
  }else{
    kappa = round(init_CP/n*N)
  }

  out = list() # An output list
  out$Kappa = matrix(data = rep(0, times = (num_iter-num_warmup)*L), nrow = num_iter-num_warmup, ncol = L)
  out$Beta = array(data = rep(x = 0, times = (num_iter-num_warmup)*p*J), dim = c((num_iter-num_warmup), p, J)) 
  out$P = array(data = rep(x = 0, times = (num_iter-num_warmup)*N*J), dim = c((num_iter-num_warmup), N, J))
  out$AUC = matrix(0, nrow = num_iter-num_warmup, ncol = J-1)

  
  #################### START of the for-loop fitting multiple change points #####################
  pb = txtProgressBar(min = 0, max = num_iter, style = 3)

  start = Sys.time()
  
  for (iter in 1:num_iter){
    
    # 1. Update samples 
    gs = GS_update(temper = 1, beta, kappa, omega, data, X, X_auc, N, L, J, m0, V0inv)
    beta = gs$beta
    kappa = gs$kappa
    omega = gs$omega
    P = gs$P
    AUC = gs$Auc
    original_kappa = gs$original_kappa
    
    
    # 2. Store the samples in the current iteration
    if ( iter > num_warmup){
      out$Beta[iter-num_warmup,,] = beta
      out$Kappa[iter-num_warmup,] = original_kappa
      out$P[iter-num_warmup,,] = P
      out$AUC[iter-num_warmup,] = AUC
    }
    
    setTxtProgressBar(pb, iter)
  }
  
  
  end = Sys.time()
  close(pb)
  
  ##################### END of the for-loop fitting multiple change points #######################
  
  runTime = end - start
  
  if (print_outputs){
    cat("\n" ,"The algorithm took ", format(runTime), "to run ", num_iter, 
        "iterations, including ", num_warmup, "burn-in iterations", "to fit ", num_CP, 
        "change points in a data set in", p, "dimensions." )
    cat("The output list has the structure", "\n", str(out))
  }
  
  return(out)
}