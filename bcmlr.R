
###################################################################################################################
# bcmlr fits single or multiple changepoints (CPs), finds an optimal tempering schedule (if fitting multiple CPs)
# and returns an output list including posterior samples of change points and regression coefficients
###################################################################################################################

bcmlr <- function(data, num_iter, num_warmup, num_temper, num_CP,  num_tune, pc_cores = detectCores()/2){
  
  # We can print a warning to users if input data is not a matrix or have NA values
  # X = as.matrix(data)  
  N = dim(X)[1] # Sample size
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
  
  
  ############################
  # SINGLE change point case #
  ############################
  if (num_CP == 1){
    
    cat("Fitting ", num_CP, " change points...... \n")
    
    # priors for beta (logistic regression coefficients)
    m0 <- rep(0, times = p) # m0[,j] is the prior mean vector 
    sd0 <- 3 # prior standard deviation 
    V0 <- sd0*diag(p) # V0 is the prior cov. matrix
    V0inv <- 1/sd0*diag(p) # V0inv is the prior precision matrix
    
    
    # Initialize working quantities
    beta = rep(0, times = p) 
    omega = rep(0, times = N) 
    P = rep(0, times = N) 
    # Initialize at the middle time point. 
    kappa = N/2 
    
    # Create an output list
    out = list()
    out$Beta = matrix(0, nrow = num_iter-num_warmup, ncol = p)
    out$P = matrix(0, nrow = num_iter-num_warmup, ncol = N)
    out$Kappa = rep(0, times = num_iter - num_warmup)
    out$AUC = rep(0, times = num_iter - num_warmup)
    
    
    ######### START of the for-loop ########
    start = Sys.time()
    pb = txtProgressBar(min = 0, max = num_iter, style = 3)
    for (iter in 1:num_iter){
      ####### Update omega ######
      eta = X %*% beta
      for (i in 1:N){
        omega[i] = rpg.devroye(num = 1, h = 1, z = eta[i]) 
      }
      
      ####### Update beta #######
      V = solve(t(X)%*%diag(omega)%*%X + V0inv)
      # Update delta based on the kappa from the previous iteration
      delta = c(rep(-0.5, times = kappa), rep(0.5, times = N-kappa)) 
      m = V%*%(t(X)%*%delta + V0inv%*%m0)
      beta = mvrnorm(n=1, mu = m, Sigma = V)
      
      ####### Update P #########
      Phi_exp <- exp(X %*% beta) 
      P = Phi_exp/(1+Phi_exp)
      
      range <- 1:(N-1)
      lr <- length(range)
      log_prob_vec_prop <- rep(0, lr) 
      for(l in 1:lr){
        log_prob_vec_prop[l] = sum(log(1-P[1:range[l]])) + sum(log(P[(range[l]+1):N]))
      }
      prob_vec <- exp(log_prob_vec_prop - logsumexp(log_prob_vec_prop)) 
      if(lr > 1){
        kappa = sample(x=range, size=1, prob=prob_vec)
      }
      
      # Update Y based on the current kappa
      Y = c(rep(0, times = kappa), rep(1, times = N-kappa))
      # Compute AUC value
      auc = calculate_roc_auc(Y, P)
      
      # Store samples after burn-in
      if(iter > num_warmup){
        out$Beta[iter-num_warmup,] = beta 
        out$P[iter-num_warmup,] = P
        out$Kappa[iter-num_warmup] = kappa
        out$AUC[iter-num_warmup] = auc
      }
      setTxtProgressBar(pb, iter)
    }
    end = Sys.time()
    close(pb)
    ######### END of the for-loop ##########
    
    runTime = end - start
  }

  ##############################
  # Multiple change point case #
  ##############################
  else{
    
    # bisection, UpdateSchedule, DEO are functions used for finding the optimal tempering schedule
    # Reference: Syed et al.(2022)   https://academic.oup.com/jrsssb/article/84/2/321/7056147
    
    
    # bisection and DEO implement Algorithm 2 in Syed et al.(2022)
    bisection <- function(f, a = 0, b = 1, tol = 1e-06){
      x = a+(b-a)/2
      while ( abs(b-a) > tol) {
        if (f(a)*f(x) < 0) {
          b <- x
        } else {
          a <- x
        }
        x <- a+(b-a)/2
      }
      return(x)
    }
    
    UpdateSchedule <- function(f, schedule.size){
      Lambda = f(1)
      T.set.tune = rep(0, times = schedule.size)
      for (k in 1:(schedule.size-1)){
        target = k/schedule.size*Lambda
        f_minus_target <- function(x) f(x) - target
        T.set.tune[k] = bisection(f = f_minus_target)
      }
      T.set.tune[schedule.size] = 1
      
      return(T.set.tune)
    }
    
    
    # DEO (deterministic even-odd) is a function that returns 
    # a list of rejection probabilities
    # and a tempering schedule adapted to the data
    DEO <- function(data, num_iter, schedule, num_CP,  pc_cores = detectCores()/2){

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
      
      # I used this function when calculating the log target in MH steps. 
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
      
      
      # https://arxiv.org/pdf/2205.04997
      # https://gregorygundersen.com/blog/2020/02/09/log-sum-exp/
      logsumexp <- function(x){
        c <- max(x)
        return(c + log(sum(exp(x - c))))
      }
      
      
      # A function that updates the MNL coefficient matrix \beta, the change point 
      # vector \kappa, and the auxiliary matrix \omega per iteration 
      GS_update <- function(temper, beta, kappa, omega, X, N, L, J, m0, V0inv){
        
        # update Y, delta 
        Y = kappa_to_Y(kappa, N, L)
        delta = temper*(Y-1/2)
        
        ########################
        # Update beta and omega#
        ########################
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
          V_j <- chol2inv(chol(t(X) %*% diag(omega[,j]) %*% X + V0inv[,,j]))
          m_j <- V_j %*% (t(X) %*% (delta[,j] + diag(omega[,j]) %*% c_j) + V0inv[,,j] %*% m0[,j])
          beta[,j] = mvrnorm(n = 1, mu = m_j, Sigma = V_j)  # update beta
          
        }
        #########################################################
        # Update kappa (and Y, which is a function of kappa) # 
        #########################################################
        
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
              log_prob_vec_prop[l] <- temper*(sum(log(P[1:range[l],j])) + sum(log(P[(range[l]+1):kappa[j+1],j+1])) - 
                                                length(1:range[l])*log(length(1:range[l])) - 
                                                length((range[l]+1):kappa[j+1])*log(length((range[l]+1):kappa[j+1])))
            }
            #prob_vec <- exp(log_prob_vec_prop)/sum(exp(log_prob_vec_prop))
            prob_vec <- exp(log_prob_vec_prop - logsumexp(log_prob_vec_prop)) 
            
          } else if(j==L){
            range <- (kappa[j-1]+1):(N-1)
            lr <- length(range)
            log_prob_vec_prop <- rep(0, lr)
            for(l in 1:lr){
              log_prob_vec_prop[l] <- temper*(sum(log(P[(kappa[j-1]+1):range[l],j])) + sum(log(P[(range[l]+1):N,j+1])) 
                                              - length((kappa[j-1]+1):range[l])*log(length((kappa[j-1]+1):range[l]))
                                              - length((range[l]+1):N)*log(length((range[l]+1):N)))
            }
            #prob_vec <- exp(log_prob_vec_prop)/sum(exp(log_prob_vec_prop))
            prob_vec <- exp(log_prob_vec_prop - logsumexp(log_prob_vec_prop)) 
          } else{
            range <- (kappa[j-1]+1):(kappa[j+1]-1)
            lr <- length(range)
            log_prob_vec_prop <- rep(0, lr) 
            for(l in 1:lr){
              log_prob_vec_prop[l] <- temper*(sum(log(P[(kappa[j-1]+1):range[l],j])) + sum(log(P[(range[l]+1):kappa[j+1],j+1])) 
                                              - length((kappa[j-1]+1):range[l])*log(length((kappa[j-1]+1):range[l])) 
                                              - length((range[l]+1):kappa[j+1])*log(length((range[l]+1):kappa[j+1])))
            }
            #prob_vec <- exp(log_prob_vec_prop)/sum(exp(log_prob_vec_prop))
            prob_vec <- exp(log_prob_vec_prop - logsumexp(log_prob_vec_prop)) 
          }
          
          if(lr > 1){
            kappa[j] <- sample(x=range, size=1, prob=prob_vec) # update kappa
          }
        }
        
        return( list(beta = beta, kappa = kappa, omega = omega) )
        
      }
      
      
      
      # A function that return a value of the logarithm of the target density kernel
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
      
      # A function that "wraps" in the Gibbs sampler update, 
      # specifies under which tempering power the Gibbs sampler is run
      GS_parallel <- function(t) {
        GS_update(T.set[t], beta[t,,], kappa[t,], omega[t,,], X, N, L, J, m0, V0inv)
      }
      
      ##################
      # A big for-loop #
      ##################
      
      N = dim(data)[1] # Sample size
      p = dim(data)[2] # data dimension
      L = num_CP # the number of change points 
      J = L+1 # the number of classes
      
      m0 <- array(data=0, dim=c(p,J-1)) # m0[,j] is the prior mean of the jth coef vec
      V0 <- array(data=0, dim=c(p,p,J-1)) # V0[,,j] is the prior cov. of the jth coef vec
      V0inv <- array(data=0, dim=c(p,p,J-1)) # V0inv[,,j] is the inverse of V0[,,j]
      sd0 <- 3 # prior standard deviation 
      for(j in 1:(J-1)){
        V0[,,j] <- sd0*diag(p)
        V0inv[,,j] <- 1/sd0*diag(p)
      }
      
      T.set = schedule
      num_temper = length(schedule)
      
      # Initialize working quantities
      beta = array(data = rep(x = 0, times = num_temper*p*J), dim = c(num_temper, p, J))
      omega = array(data = rep(x = 0, times = num_temper*N*J), dim = c(num_temper, N, J))
      kappa = matrix(data = 0, nrow = num_temper, ncol = L)
      for ( r in 1:num_temper){
        for (l in 1:L){
          kappa[r, l] = round(N/J) *l
        }
      }
      rej = rep(0, times = num_temper - 1)
      
      
      indices = 1:(num_temper-1) 
      Ind.even = indices[indices%%2 == 0] # even indices 
      Ind.odd = indices[indices%%2 != 0] # odd indices
      
      for (iter in 1:num_iter){
        
        # 1. Update samples at all tempering powers
        # parallel computing:
        gs = mclapply(X = 1:num_temper, FUN = GS_parallel, mc.cores = pc_cores)
        # store the results into the working/computation quantities: 
        for (t in 1:num_temper){
          beta[t,,] = gs[[t]]$beta 
          kappa[t,] = gs[[t]]$kappa
          omega[t,,] = gs[[t]]$omega
        }
        
        
        # 2. Determine even/odd indices based on the current iteration 
        if ( iter %%2 == 0){
          Ind.iter = Ind.even 
        }else{
          Ind.iter = Ind.odd
        }
        
        # 3. Do Metropolis-Hastings checks and potential swaps on all pairs of samples in the current E/O index set 
        
        # ***Parallel computing (if all pairs are completely independent)***
        # Under the DEO scheme, all pairs are indeed independent per iteration. 
        # mh = mclapply(X = Ind.iter, FUN = MH_parallel, mc.cores = num_cores)
        # pair_swap = mh$pair_swap
        
        # ***Serial computing***
        for (id in Ind.iter){
          t1 = id
          t2 = id + 1
          
          p12 = logtarget(T.set[t1], beta[t2,,], kappa[t2,], omega[t2,,], X, N, L, J, m0, V0inv)
          p21 = logtarget(T.set[t2], beta[t1,,], kappa[t1,], omega[t1,,], X, N, L, J, m0, V0inv)
          p11 = logtarget(T.set[t1], beta[t1,,], kappa[t1,], omega[t1,,], X, N, L, J, m0, V0inv)
          p22 = logtarget(T.set[t2], beta[t2,,], kappa[t2,], omega[t2,,], X, N, L, J, m0, V0inv)
          
          logA = p12 + p21 - p11 - p22
          u = log(runif(n = 1))
          rho = min(0,logA)
          
          rej[id] = rej[id] + 1 - exp(rho) 
          
          if ( u < rho ){
            
            beta_t1 = beta[t1,,]
            kappa_t1 = kappa[t1,]
            omega_t1 = omega[t1,,]
            
            beta[t1,,] = beta[t2,,]
            kappa[t1,] = kappa[t2,]
            omega[t1,,] = omega[t2,,]
            
            beta[t2,,] = beta_t1
            kappa[t2,] = kappa_t1
            omega[t2,,] = omega_t1
          }
        }
        
      }
      
      ave_rej = rej/num_iter
      
      return(ave_rej)
    }
    
    #################### Finding the optimal tempering schedule ###########################
    
    cat("Finding the optimal tempering schedule......", "\n")
    
    # Initial tempering schedule with a user-specified length
    T.set = seq(from = 0.0001, to = 1, length = num_temper)
    # Tuning iterations used to adjust the tempering schedule
    num_tune =  num_tune
    maxRound = log2(num_tune)
    
    # Iterations (scans) per round, increased exponentially
    n = 1  
    # a vector of rejection probabilities
    rej = rep(0, times = num_temper - 1)
    
    start = Sys.time()
    pb = txtProgressBar(min = 0, max = maxRound, style = 3)
    for (rd in 1:maxRound){
      
      # 1. Rejection probabilities
      rej = DEO(data = X, num_iter = n, num_CP = num_CP, schedule = T.set)
      
      # 2. Update the Communication Barrier, Lambda.fun: 
      # For each tempering power, compute Lambda as in Eq.(32) in Syed et al. (2021)
      # After the for-loop, compute a monotone increasing interpolation. 
      Lambda.vec = rep(0, times = num_temper)
      for ( k in 1:(num_temper-1)){
        Lambda.vec[k+1] = sum(rej[1:k])
      }
      Lambda.fun = splinefun(T.set, Lambda.vec, method = "monoH.FC")
      
      # 3. UpdateSchedule
      T.set = UpdateSchedule(f = Lambda.fun, schedule.size = num_temper)
      
      n = 2*n # rounds use an exponentially increasing number of scans
      setTxtProgressBar(pb, rd)
    }
    close(pb)
    end = Sys.time()
    
    # Time used to find the optimal tempering schedule
    runTime = end - start
    
    # Optimal schedule size (number of tempering powers)
    optim.size = 2*Lambda.fun(1)
    # Optimal tempering schedule
    T.set = UpdateSchedule(f = Lambda.fun, schedule.size = optim.size)
    # Update num_temper to the optimal size = length(T.set)
    num_temper = optim.size
    
    cat("It took ", format(runTime), "to find the optimal tempering schedule.")

    
    #########################################################################################################
    # Utility functions and functions for Gibbs sampling, Metropolis-Hastings check, and parallel computing # 
    #########################################################################################################
    
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
  
    
    # A function that updates the MNL coefficient matrix \beta, the change point 
    # vector \kappa, and the auxiliary matrix \omega per iteration 
    GS_update <- function(temper, beta, kappa, omega, X, N, L, J, m0, V0inv){
      
      # update Y, delta 
      Y = kappa_to_Y(kappa, N, L)
      delta = temper*(Y-1/2)
      
      ########################
      # Update beta and omega#
      ########################
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
        V_j <- chol2inv(chol(t(X) %*% diag(omega[,j]) %*% X + V0inv[,,j]))
        m_j <- V_j %*% (t(X) %*% (delta[,j] + diag(omega[,j]) %*% c_j) + V0inv[,,j] %*% m0[,j])
        beta[,j] = mvrnorm(n = 1, mu = m_j, Sigma = V_j)  # update beta
        
      }
      #########################################################
      # Update kappa (and Y, which is a function of kappa) # 
      #########################################################
      
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
            log_prob_vec_prop[l] <- temper*(sum(log(P[1:range[l],j])) + sum(log(P[(range[l]+1):kappa[j+1],j+1])) - 
                                              length(1:range[l])*log(length(1:range[l])) - 
                                              length((range[l]+1):kappa[j+1])*log(length((range[l]+1):kappa[j+1])))
          }
          #prob_vec <- exp(log_prob_vec_prop)/sum(exp(log_prob_vec_prop))
          prob_vec <- exp(log_prob_vec_prop - logsumexp(log_prob_vec_prop)) 
          
        } else if(j==L){
          range <- (kappa[j-1]+1):(N-1)
          lr <- length(range)
          log_prob_vec_prop <- rep(0, lr)
          for(l in 1:lr){
            log_prob_vec_prop[l] <- temper*(sum(log(P[(kappa[j-1]+1):range[l],j])) + sum(log(P[(range[l]+1):N,j+1])) 
                                            - length((kappa[j-1]+1):range[l])*log(length((kappa[j-1]+1):range[l]))
                                            - length((range[l]+1):N)*log(length((range[l]+1):N)))
          }
          #prob_vec <- exp(log_prob_vec_prop)/sum(exp(log_prob_vec_prop))
          prob_vec <- exp(log_prob_vec_prop - logsumexp(log_prob_vec_prop)) 
        } else{
          range <- (kappa[j-1]+1):(kappa[j+1]-1)
          lr <- length(range)
          log_prob_vec_prop <- rep(0, lr) 
          for(l in 1:lr){
            log_prob_vec_prop[l] <- temper*(sum(log(P[(kappa[j-1]+1):range[l],j])) + sum(log(P[(range[l]+1):kappa[j+1],j+1])) 
                                            - length((kappa[j-1]+1):range[l])*log(length((kappa[j-1]+1):range[l])) 
                                            - length((range[l]+1):kappa[j+1])*log(length((range[l]+1):kappa[j+1])))
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
      Y <- kappa_to_Y(kappa, N, L)
      auc_values = rep(0, times = J-1)
      for (j in 1:(J-1)){
        if(j == 1){
          rg = 1:(kappa[j+1]-1)
          ratios = P[rg,j+1]/(P[rg,j] + P[rg,j+1])
          auc_values[j] = calculate_roc_auc(Y[rg, j + 1], ratios)
        }else if (j == J-1){
          rg = (kappa[j-1]+1):(N-1)
          ratios = P[rg,j+1]/(P[rg,j] + P[rg, j+1])
          auc_values[j] = calculate_roc_auc(Y[rg, j + 1], ratios)
        }else{
          rg = (kappa[j-1]+1):(kappa[j+1]-1)
          ratios = P[rg,j+1]/(P[rg,j] + P[rg, j+1])
          auc_values[j] = calculate_roc_auc(Y[rg, j + 1], ratios)
        }
      }
      
      if (temper == 1){
        return( list(beta = beta, kappa = kappa, omega = omega, P = P, Auc = auc_values) )
      }else{
        return( list(beta = beta, kappa = kappa, omega = omega) )
      }
      
    }
    
    # A function that return a value of the logarithm of the target density kernel
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
    
    
    # A function that "wraps" in the Gibbs sampler update, 
    # specifies under which tempering power the Gibbs sampler is run
    GS_parallel <- function(t) {
      GS_update(T.set[t], beta[t,,], kappa[t,], omega[t,,], X, N, L, J, m0, V0inv)
    }
    
    
    ##################################### Fitting multiple change points ####################################
    
    cat("\n", "Fitting ", num_CP," change points with the tempering schedule ", T.set, "...... \n")
    
    L = num_CP # the number of change points 
    J = L+1 # the number of classes
    
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
    beta = array(data = rep(x = 0, times = num_temper*p*J), dim = c(num_temper, p, J))
    omega = array(data = rep(x = 0, times = num_temper*N*J), dim = c(num_temper, N, J))
    P = matrix(0, nrow = N, ncol = J)
    Auc = rep(0, times = J-1)
    kappa = matrix(data = 0, nrow = num_temper, ncol = L)
    for ( r in 1:num_temper){
      for (l in 1:L){
        kappa[r, l] = round(N/J) *l
      }
    }
    rej = rep(0, times = num_temper - 1) 
    
    
    out = list() # An output list
    out$Kappa = array(data = rep(x = 0, times = (num_iter-num_warmup)*num_temper*L), dim = c((num_iter-num_warmup), num_temper, L))
    # The following three lists only store samples under tempering power 1. 
    out$Beta = array(data = rep(x = 0, times = (num_iter-num_warmup)*p*J), dim = c((num_iter-num_warmup), p, J)) 
    out$P = array(data = rep(x = 0, times = (num_iter-num_warmup)*N*J), dim = c((num_iter-num_warmup), N, J))
    out$AUC = matrix(0, nrow = num_iter-num_warmup, ncol = J-1)
    
    pair_swap = rep(0, times = num_temper-1)
    
    indices = 1:(num_temper-1) 
    Ind.even = indices[indices%%2 == 0] # even indices 
    Ind.odd = indices[indices%%2 != 0] # odd indices
    
    #################### START of the for-loop fitting multiple change points #####################
    pb = txtProgressBar(min = 0, max = num_iter, style = 3)
    start = Sys.time()
    for (iter in 1:num_iter){
      
      # 1. Update samples at all tempering powers
      # parallel computing:
      gs = mclapply(X = 1:num_temper, FUN = GS_parallel, mc.cores = pc_cores)
      
      # store the results into the working/computation quantities: 
      for (t in 1:num_temper){
        kappa[t,] =  gs[[t]]$kappa
        beta[t,,] =  gs[[t]]$beta 
        omega[t,,] = gs[[t]]$omega
      }
      P = gs[[num_temper]]$P
      Auc = gs[[num_temper]]$Auc
      
      
      
      # 2. Determine even/odd indices based on the current iteration 
      if ( iter %%2 == 0){
        Ind.iter = Ind.even 
      }else{
        Ind.iter = Ind.odd
      }
      
      # 3. Do Metropolis-Hastings checks and potential swaps on all pairs of samples in the current E/O index set 
      
      # ***Parallel computing (if all pairs are completely independent)***
      # Under the DEO scheme, all pairs are indeed independent per iteration. 
      # mh = mclapply(X = Ind.iter, FUN = MH_parallel, mc.cores = num_cores)
      # pair_swap = mh$pair_swap
      
      # ***Serial computing***
      for (id in Ind.iter){
        t1 = id
        t2 = id + 1
        
        p12 = logtarget(T.set[t1], beta[t2,,], kappa[t2,], omega[t2,,], X, N, L, J, m0, V0inv)
        p21 = logtarget(T.set[t2], beta[t1,,], kappa[t1,], omega[t1,,], X, N, L, J, m0, V0inv)
        p11 = logtarget(T.set[t1], beta[t1,,], kappa[t1,], omega[t1,,], X, N, L, J, m0, V0inv)
        p22 = logtarget(T.set[t2], beta[t2,,], kappa[t2,], omega[t2,,], X, N, L, J, m0, V0inv)
        
        logA = p12 + p21 - p11 - p22
        u = log(runif(n = 1))
        rho = min(0,logA)
        
        rej[id] = rej[id] + 1 - exp(rho) 
        
        if ( u < rho ){
          
          beta_t1 = beta[t1,,]
          kappa_t1 = kappa[t1,]
          omega_t1 = omega[t1,,]
          
          beta[t1,,] = beta[t2,,]
          kappa[t1,] = kappa[t2,]
          omega[t1,,] = omega[t2,,]
          
          beta[t2,,] = beta_t1
          kappa[t2,] = kappa_t1
          omega[t2,,] = omega_t1
          
          pair_swap[id] = pair_swap[id] + 1
        }
      }
      
      
      # 4. Store the samples in the current iteration
      if ( iter > num_warmup){
        out$Beta[iter-num_warmup,,] = beta[num_temper,,]
        out$Kappa[iter-num_warmup,,] = kappa
        out$P[iter-num_warmup,,] = P
        out$AUC[iter-num_warmup,] = Auc
      }
      
      setTxtProgressBar(pb, iter)
    }
    
    end = Sys.time()
    close(pb)
    
    ##################### END of the for-loop fitting multiple change points #######################
    
    
    
    # Store the average MH acceptance rates, posterior mean rejection probabilities, 
    # and the optimal tempering schedule used for this algorithm in the output list.
    out$pair_acc = pair_swap/(num_iter/2)
    out$Rej = rej/num_iter
    out$schedule = T.set
    
    runTime = end - start
  }
  
  cat("\n" ,"The algorithm took ", format(runTime), "to run ", num_iter, 
      "iterations, including ", num_warmup, "burn-in iterations", "to fit ", num_CP, 
      "change points in a data set in", p, "dimensions.", "\n" )
  cat("The output list has the structure", str(out))
  
  return(out)
}
