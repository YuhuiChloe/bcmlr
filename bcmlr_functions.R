#############
# Libraries #
#############

library(BayesLogit) # for rpg
library(MASS) # for mvrnorm
library(PGdensityDraft) # for PG density
# library(mvtnorm) # for dmvnorm

library(parallel) # for parallel computing using mclapply()

library(ggplot2)  # for pair plots
library(GGally)

library(knitr) # for printing beta matrix using kable()



####################################################################################
# A function that runs a for-loop with Gibbs samplers and Metropolis-Hastings checks 
# and prints results #
####################################################################################
# bcmlr <- function(data, T.set, P, beta, kappa, omega, X, N, L, J, m0, V0inv) 

bcmlr <- function(data,  num_CP = 3, num_temper = 30, num_iter = 10000,  num_warmup = 5000,   num_cores = detectCores()/2){
  
  
  #####################
  # Utility functions # 
  #####################
  
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
  # A function that converts a response matrix with the change point structure 
  # to a change point vector
  Y_to_kappa <- function(Y, L){
    # L <- dim(Y)[2] - 1
    kappa <- rep(NA, L)
    for(l in 1:L){
      kappa[l] <- max(which(Y[,l] == 1))
    }
    return(kappa)
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
  
  
  ##############################################################
  # Functions for Gibbs sampling and Metropolis-Hastings check #
  ##############################################################
  
  # A function that updates the MNL coefficient matrix \beta, the change point 
  # vector \kappa, and the auxilliary matrix \omega per iteration 
  GS_update <- function(temper, P, beta, kappa, omega, X, N, L, J, m0, V0inv){
    
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
    
    return( list(beta = beta, kappa = kappa, omega = omega, P = P) )
    
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
  
  ####################################
  # Functions for parallel computing #
  ####################################
  
  # A function that "wraps" in the Gibbs sampler update, 
  # specifies under which tempering power the Gibbs sampler is run
  GS_parallel <- function(t) {
    GS_update(T.set[t], P[t,,], beta[t,,], kappa[t,], omega[t,,], X, N, L, J, m0, V0inv)
  }
  
  
  
  
  

  
  
  # X = as.matrix(data)  
  N = dim(X)[1] # Sample size
  p = dim(X)[2] # data dimension
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
  
  T.set = seq(.05^(1/4), 1, length.out = num_temper)^4 # The default set of tempering powers
  
  # initialize working quantities
  beta = array(data = rep(x = 0, times = num_temper*p*J), dim = c(num_temper, p, J))
  omega = array(data = rep(x = 0, times = num_temper*N*J), dim = c(num_temper, N, J))
  P = array(data = rep(x = 0, times = num_temper*N*J), dim = c(num_temper, N, J))
  kappa = matrix(data = 0, nrow = num_temper, ncol = L)
  for ( r in 1:num_temper){
    for (l in 1:L){
      kappa[r, l] = round(N/J) *l
    }
  }
  
  out = list() # An output list
  out$Beta = array(data = rep(x = 0, times = (num_iter-num_warmup)*p*J), dim = c(num_iter-num_warmup, p, J))
  out$Kappa = array(data = rep(x = 0, times = (num_iter-num_warmup)*num_temper*L), dim = c(num_iter-num_warmup, num_temper, L))
  out$P = array(data = rep(x = 0, times = (num_iter-num_warmup)*N*J), dim = c(num_iter-num_warmup, N, J))
  
  MH_accept = 0
  pair_swap = rep(0, times = num_temper-1)
  
  indices = 1:(num_temper-1) 
  Ind.even = indices[indices%%2 == 0] # even indices 
  Ind.odd = indices[indices%%2 != 0] # odd indices
  
  pb = txtProgressBar(min = 0, max = num_iter, style = 3)
  start = Sys.time()
  for (iter in 1:num_iter){
    
    # 1. Update samples at all tempering powers
    # parallel computing:
    gs = mclapply(X = 1:num_temper, FUN = GS_parallel, mc.cores = num_cores)

    # store the results into the working/computation quantities: 
    for (t in 1:num_temper){
      beta[t,,] = gs[[t]]$beta 
      kappa[t,] = gs[[t]]$kappa
      omega[t,,] = gs[[t]]$omega
      P[t,,] = gs[[t]]$P
    }
    
    # 2. Determine even/odd indices based on the current iteration 
    if ( iter %%2 ==0){
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
        MH_accept = MH_accept + 1
      }
    }
    
    
    # 4. Store the samples in the current iteration
    if ( iter > num_warmup){
      out$Beta[iter-num_warmup,,] = beta[num_temper,,]
      out$P[iter-num_warmup,,] = P[num_temper,,]
      out$Kappa[iter-num_warmup,,] = kappa
    }
    
    # # 5. Print progress
    # if (iter%%1000 == 0){
    #   print(c(iter, MH_accept/(iter*(num_temper-1)/2)))
    #   #print(pair_swap/(iter/2))
    # }
    
    setTxtProgressBar(pb, iter)
  }
  
  end = Sys.time()
  close(pb)
  
  runTime = end - start
  paste("The algorithm took ", format(runTime), "to run ", num_iter, 
        "iterations, including ", num_warmup, "burn-in iterations.")
  
  return(out)
}
