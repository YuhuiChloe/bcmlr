# bcmlr fits single or multiple changepoints (CPs), finds an optimal tempering schedule (if fitting multiple CPs)
# and returns an output list including posterior samples of change points and regression coefficients
# Input data needs to be in a matrix form.
# num_CP: the number of change points to fit
# num_tune:  determine the number of tuning iterations needed to adaptively find an optimal tempering schedule
# num_temper: the number of tempering powers = the tempering schedule length
# num_iter: number of iterations
# num_warmup: number of burn-in 
# pc_cores: the number of cores used for parallel computing

bcmlr_PT <- function(data, num_iter = 5000, num_warmup = 2500, num_tune = 5000, num_temper = 40, num_CP,  pc_cores = detectCores(), print_outputs = TRUE){
  
  X = as.matrix(data)
  N = dim(X)[1] # Sample size
  p = dim(X)[2] # data dimension
  #### Functions used for both cases (single CP and multiple CPs) ####
  # https://arxiv.org/pdf/2205.04997
  # https://gregorygundersen.com/blog/2020/02/09/log-sum-exp/
  # A function that mitigate computational issues in 
  # massive division or multiplication of extreme values. 
  logsumexp <- function(x){
    c <- max(x)
    return(c + log(sum(exp(x - c))))
  }
  
  ##### SINGLE change point case: no PT needed ######
  if (num_CP == 1){
    cat("Fitting ", num_CP, " change points...... \n")
    # priors for beta (logistic regression coefficients)
    m0 <- rep(0, times = p) # m0[,j] is the prior mean vector 
    sd0 <- 3 # prior variance
    V0 <- sd0*diag(p) # V0 is the prior cov. matrix
    V0inv <- 1/sd0*diag(p) # V0inv is the prior precision matrix
    # Initialize working quantities
    beta = rep(0, times = p) 
    omega = rep(0, times = N) 
    P = rep(0, times = N) 
    kappa = round(N/2) #init_CP
    # Create an output list
    out = list()
    out$Beta = matrix(0, nrow = num_iter-num_warmup, ncol = p)
    out$P = matrix(0, nrow = num_iter-num_warmup, ncol = N)
    out$Kappa = rep(0, times = num_iter - num_warmup)
    out$AUC = rep(0, times = num_iter - num_warmup)
    ######## START of the for-loop
    start = Sys.time()
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
      # Update Y based on the current kappa
      Y = c(rep(0, times = kappa), rep(1, times = N-kappa))
      # Store samples after burn-in
      if(iter > num_warmup){
        out$Beta[iter-num_warmup,] = beta 
        out$P[iter-num_warmup,] = P
        out$Kappa[iter-num_warmup] = kappa
      }
      setTxtProgressBar(pb, iter)
    }
    end = Sys.time()
    close(pb) 
    # Use posterior modes as changepoint estimates
    cp = table(out$Kappa)
    out$Kappa_mode = as.integer(names(cp)[which.max(cp)]) # posterior mode 
    ######### END of the for-loop 
    if (print_outputs){
      runTime = end - start
      cat("The algorithm took ", format(runTime), "to run ", num_iter, "iterations with", 
          num_warmup, "burn-in iterations to fit",  num_CP, 
          "change point in a data set in", p, "dimensions.")
      cat("The output list has the structure", "\n", str(out))
    }
  } # end of the single CP case
  
  ####### Multiple change point case: use parallel tempering #######
  else{
    # Reference: Syed et al.(2022)   https://academic.oup.com/jrsssb/article/84/2/321/7056147
    # bisection, UpdateSchedule, DEO are functions used for finding the optimal tempering schedule
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
    # Algo 2: UpdateSchedule in Syed et al.(2022)
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
    # Algo 1: DEO in Syed et al.(2022)
    # DEO (deterministic even-odd) is a function 
    # that returns a list of rejection probabilities and a tempering schedule adapted to the data
    DEO <- function(data, num_iter, num_warmup, schedule, num_CP, num_cores = pc_cores, optimal_schedule = FALSE){
      # kappa_to_Y is a function that converts a change point vector to a response matrix
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
      # Use kappa_to_alpha when calculating the log target in MH steps. 
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
        # 0. update Y, delta 
        Y = kappa_to_Y(kappa, N, L)
        delta = temper*(Y-1/2)
        # 1. Update beta and omega
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
        # 2. Update kappa (and Y, which is a function of kappa)
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
        return( list(beta = beta, kappa = kappa, omega = omega, P = P))
      }
      # logtarget is a function that return a value of the logarithm of the target density kernel
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
      # Gibbs sampler for-loop starts below
      N = dim(data)[1] # Sample size
      p = dim(data)[2] # data dimension
      L = num_CP # the number of change points 
      J = L+1 # the number of classes
      m0 <- array(data=0, dim=c(p,J-1)) # m0[,j] is the prior mean of the jth coef vec
      V0 <- array(data=0, dim=c(p,p,J-1)) # V0[,,j] is the prior cov. of the jth coef vec
      V0inv <- array(data=0, dim=c(p,p,J-1)) # V0inv[,,j] is the inverse of V0[,,j]
      sd0 <- 3 # prior variance
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
          kappa[r, l] = round(N/J) *l # even initialization
        }
      }
      out = list()
      # The following lists store samples under tempering power 1.
      out$Kappa = matrix(data = rep(0, times = (num_iter-num_warmup)*L), nrow = num_iter-num_warmup, ncol = L)
      out$Beta = array(data = rep(x = 0, times = (num_iter-num_warmup)*p*J), dim = c((num_iter-num_warmup), p, J))
      out$P = array(data = rep(x = 0, times = (num_iter-num_warmup)*N*J), dim = c((num_iter-num_warmup), N, J))
      # a vector to store the frequency of a pair of samples being swapped at each index
      pair_swap = rep(0, times = num_temper-1)
      rej = rep(0, times = num_temper - 1)
      indices = 1:(num_temper-1) 
      Ind.even = indices[indices%%2 == 0] # even indices 
      Ind.odd = indices[indices%%2 != 0] # odd indices
      if (optimal_schedule){
        cat("\n", "Fitting ", num_CP ," change points with ", num_temper ,"tempering powers ", T.set, "...... \n")
        pb = txtProgressBar(min = 0, max = num_iter, style = 3)
        start = Sys.time() # start time
      }
      for (iter in 1:num_iter){
        # 1. Update samples at all tempering powers
        # Pre-run GS_parallel once in the parent process to initialize
        # packages/compiled code and avoid intermittent fork errors in mclapply()
        invisible(GS_parallel(1))
        # parallel computing:
        gs = mclapply(X = 1:num_temper, FUN = GS_parallel, mc.cores = detectCores())
        # store the results into the working/computation quantities: 
        for (t in 1:num_temper){
          beta[t,,] = gs[[t]]$beta 
          kappa[t,] = gs[[t]]$kappa
          omega[t,,] = gs[[t]]$omega
        }
        P = gs[[num_temper]]$P
        # 2. Determine even/odd indices based on the current iteration 
        if ( iter %%2 == 0){
          Ind.iter = Ind.even 
        }else{
          Ind.iter = Ind.odd
        }
        # 3. Do Metropolis-Hastings checks and potential swaps on all pairs of samples in the current E/O index set 
        # ***Serial computing***
        for (id in Ind.iter){
          # two temperature/tempeirng powers under which the samples are checked
          t1 = id
          t2 = id + 1
          # compute density values
          p12 = logtarget(T.set[t1], beta[t2,,], kappa[t2,], omega[t2,,], X, N, L, J, m0, V0inv)
          p21 = logtarget(T.set[t2], beta[t1,,], kappa[t1,], omega[t1,,], X, N, L, J, m0, V0inv)
          p11 = logtarget(T.set[t1], beta[t1,,], kappa[t1,], omega[t1,,], X, N, L, J, m0, V0inv)
          p22 = logtarget(T.set[t2], beta[t2,,], kappa[t2,], omega[t2,,], X, N, L, J, m0, V0inv)
          logA = p12 + p21 - p11 - p22 # log of the MH ratio
          u = log(runif(n = 1)) # generate uniform(0,1)
          rho = min(0,logA) # log of the acceptance probability
          rej[id] = rej[id] + 1 - exp(rho)  # sum rejection probabilities
          if ( u < rho ){
            #store
            beta_t1 = beta[t1,,]
            kappa_t1 = kappa[t1,]
            omega_t1 = omega[t1,,]
            #swap
            beta[t1,,] = beta[t2,,]
            kappa[t1,] = kappa[t2,]
            omega[t1,,] = omega[t2,,]
            #update
            beta[t2,,] = beta_t1
            kappa[t2,] = kappa_t1
            omega[t2,,] = omega_t1
            # record the frequency of a pair of samples being swapped at index id
            pair_swap[id] = pair_swap[id] + 1
          }
        }
        # 4. Store the samples in the current iteration
        if ( iter > num_warmup){
          out$Beta[iter-num_warmup,,] = beta[num_temper,,]
          out$Kappa[iter-num_warmup,] = kappa[num_temper,]
          out$P[iter-num_warmup,,] = P
        }
        if (optimal_schedule){
          setTxtProgressBar(pb, iter)
        }
      }
      if(optimal_schedule){
        close(pb)
        end = Sys.time() # end time 
      }
      # Use posterior modes as changepoint estimates
      cps = c()
      for (l in 1:L){
        cp = table(out$Kappa[,l])
        kappa = as.integer(names(cp)[which.max(cp)])
        cps = c(cps, kappa)
      }
      out$Kappa_mode = cps # posterior modes
      out$ave_rej = rej/num_iter # take the average rejection probabilities
      out$pair_acc = pair_swap/(num_iter/2) # average MH acceptance rates
      out$schedule = T.set #the optimal tempering schedule used for this algorithm 
      return(out)
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
      rej = DEO(data = X, num_iter = n, num_warmup = 0, schedule = T.set, num_CP = num_CP, num_cores = pc_cores, optimal_schedule = FALSE)$ave_rej
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
    runTime = end - start    # Time used to find the optimal tempering schedule
    # Optimal schedule size (number of tempering powers)
    optim.size = round(2*Lambda.fun(1))  
    # Optimal tempering schedule
    T.set = UpdateSchedule(f = Lambda.fun, schedule.size = optim.size)
    # Update num_temper to the optimal size = length(T.set)
    num_temper = optim.size
    cat("It took ", format(runTime), "to find the optimal tempering schedule with length ", num_temper)
    ## Fit bcmlr with the optimal tempering schedule ##
    out = DEO(data = X, num_iter = num_iter, num_warmup = num_warmup, num_CP = num_CP, schedule = T.set, optimal_schedule = TRUE)
  } # end of the multiple CP case
  return(out)
}
