#Held-out strategy:
  #- Hold out every samples at indices multiple of 3 (thinning = 3) when fitting bclr.  
  #- Use the held-out data to compute AUC. 

bclr <- function(data,  num_iter = 10000, num_warmup = 5000, init_CP = NA, thinning = 1, print_outputs = FALSE){
  
  # a function to center a matrix by a vector of means
  center <- function(x, m){
    for (i in 1:ncol(x)){
      x[,i] = x[,i] - m[i]
    }
    return(x)
  }
  
  n = dim(data)[1] # length of the entire series
  X = as.matrix(data)
  
  if (length(thinning) > 2 || thinning < 0 ){
      print("thinning must be a positive integer")
  }

  if (thinning != 1){
    X_auc = X[seq_len(n) %% thinning == 0, ] # rows of X selected to calculate AUC
    X = X[seq_len(n) %% thinning != 0, ] # rows of X selected to fit bclr 
    m = colMeans(X) # column means in the subseries used to fit the model
    X = center(X, m) # center the subseries used to fit the model with its column means
    X_auc = center(X_auc, m) # center the held-out samples with the same column means 
  }else{
    m = colMeans(X) 
    X = center(X, m)
    X_auc = X
  }
  
  N = dim(X)[1] # Sample size for fitting the model
  p = dim(X)[2] # data dimension
  
  
  ##############################################################
  # Functions used for both cases (single CP and multiple CPs) #
  ##############################################################
  
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
  
  
  
  ##### SINGLE change point case ######
  #cat("Fitting ", num_CP, " change points...... \n")
  
  # priors for beta (logistic regression coefficients)
  m0 <- rep(0, times = p) # m0[,j] is the prior mean vector 
  sd0 <- 3 # prior standard deviation 
  V0 <- sd0*diag(p) # V0 is the prior cov. matrix
  V0inv <- 1/sd0*diag(p) # V0inv is the prior precision matrix
  
  
  # Initialize working quantities
  beta = rep(0, times = p) 
  omega = rep(0, times = N) 
  P = rep(0, times = N) 

  # Initialize a changepoint
  if (is.na(init_CP)){
    # if there's no user-specified initial CPs, initialize at the middle point
    kappa = round(N/2) 
  }else{
    if (thinning == 1){
      kappa = init_CP
    }else{
      # time points selected to fit the model
      thin_idx = seq_len(n)[seq_len(n) %% thinning != 0]
      if (init_CP %in% thin_idx){
        # if the initial CP is selected to fit the model, use it as is. 
        kappa = which(thin_idx == init_CP) #init_CP
      }else{
        # if the initial CP is not selected to fit the model, 
        # pick the point closest before the initial CP in the selected points. 
        left_idx = thin_idx[thin_idx < init_CP]
        if (length(left_idx) > 0){
          kappa = which(thin_idx == max(left_idx)) #max(left_idx)
        }else{
          # if there's not points to the left of the initial CP, 
          # (print a warning and) choose the middle point as initial.
          # print("no points before the initial CP")
          kappa = round(N/2)
        }
      }
    }
  }
  print(c(init_CP, kappa))
  
  # Create an output list
  out = list()
  out$Beta = matrix(0, nrow = num_iter-num_warmup, ncol = p)
  out$P = matrix(0, nrow = num_iter-num_warmup, ncol = N)
  out$Kappa = rep(0, times = num_iter - num_warmup)
  out$AUC = rep(0, times = num_iter - num_warmup)
  
  
  ######## START of the for-loop
  start_time = Sys.time()
  pb = txtProgressBar(min = 0, max = num_iter, style = 3)
  for (iter in 1:num_iter){
    ####### Update omega 
    eta = X %*% beta
    for (i in 1:N){
      omega[i] = rpg.devroye(num = 1, h = 1, z = eta[i]) 
    }
    
    ####### Update beta 
    # V = solve(t(X)%*%diag(omega)%*%X + V0inv)
    # delta = c(rep(-0.5, times = kappa), rep(0.5, times = N-kappa)) 
    # m = V%*%(t(X)%*%delta + V0inv%*%m0)
    # beta = mvrnorm(n=1, mu = m, Sigma = V)
    delta = c(rep(-0.5, times = kappa), rep(0.5, times = N-kappa)) 
    Z <- matrix(rnorm(p), p, 1) # A p-vector of standard normals
    U <- chol(crossprod(x = X, y = diag(omega) %*% X) + V0inv) # Upper Cholesky factor of the inverse of V_j 
    beta <- chol2inv(U) %*% (t(X) %*%delta + V0inv %*% m0) + backsolve(U, Z)
    ####### Update P 
    Phi_exp <- exp(X %*% beta) 
    P = Phi_exp/(1+Phi_exp)
    
    range <- 1:(N-1)
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
    
    ##### Compute AUC ####
    all_idx = 1:n
    if (thinning == 1){
      thin_idx = all_idx
    }else{
      thin_idx = all_idx[all_idx%%thinning != 0]
    }
    original_kappa = thin_idx[kappa]
    # Update Y based on the kappa in the entire series
    Y = c(rep(0, times = original_kappa), rep(1, times = n - original_kappa))
    Y_auc = Y[seq_len(length(Y)) %% thinning == 0]
    # Update the matrix of probabilities of success based on the fitted coefficients
    Phi_exp_auc <- exp(X_auc %*% beta) 
    P_auc = Phi_exp_auc/(1+Phi_exp_auc)
    # Compute AUC value
    auc = calculate_roc_auc(Y_auc, P_auc)
    
    # Store samples after burn-in
    if(iter > num_warmup){
      out$Beta[iter-num_warmup,] = beta 
      out$P[iter-num_warmup,] = P
      out$Kappa[iter-num_warmup] = original_kappa
      out$AUC[iter-num_warmup] = auc
    }
    setTxtProgressBar(pb, iter)
  }
  end_time  = Sys.time()
  close(pb) 

  ######### END of the for-loop 
  if (print_outputs){
    runTime = end_time  - start_time 
    cat("The algorithm took ", format(runTime), "to run ", num_iter, "iterations with", 
        num_warmup, "burn-in iterations to fit",  num_CP, 
        "change point in a data set in", p, "dimensions.")
    cat("The output list has the structure", "\n", str(out))
  }
  
  
  return(out)
}
