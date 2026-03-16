###################### A function to simulate X and psi ###########################
# A function that takes no inputs,
# simulates multivariate Gaussian X with TWO changes in covariances
# and returns the 2-degree polynomial embedding X
sim_mvn_cicm <- function(N, kappa0){
  v = 1
  s = 0.7
  p = 4
  Sig1 = matrix(data = c(v, s, 0, s,
                         s, v, 0, 0,
                         0, 0, v, 0,
                         s, 0, 0, v), nrow = p, ncol = p, byrow = TRUE)
  X1 = mvrnorm(n = kappa0[1], mu = rep(0,p), Sigma = Sig1)
  X2 = mvrnorm(n = kappa0[2]- kappa0[1], mu = rep(1,p), Sigma = Sig1)
  X3 = mvrnorm(n = N - kappa0[2], mu = rep(1,p), Sigma = diag(p))
  X = rbind(X1, X2, X3)
  x1 = X[,1]
  x2 = X[,2]
  x3 = X[,3]
  x4 = X[,4]
  x1_sq = X[,1]^2
  x2_sq = X[,2]^2
  x3_sq = X[,3]^2
  x4_sq = X[,4]^2
  x1x2 = X[,1]*X[,2]
  x1x3 = X[,1]*X[,3]
  x1x4 = X[,1]*X[,4]
  x2x3 = X[,2]*X[,3]
  x2x4 = X[,2]*X[,4]
  x3x4 = X[,3]*X[,4]
  # 2-degree polynomial feature embedding of X = (x1, x2, x3, x4)
  psi = cbind(x1, x2, x3, x4,
              x1_sq, x2_sq, x3_sq, x4_sq,
              x1x2, x1x3, x1x4, x2x3, x2x4, x3x4)
  data = list()
  data$X = X
  data$psi = psi
  return(data)
}


sim_mvn_cicm <- function(N, kappa0){
  v = 1
  s = 0.9
  p = 8
  Sig1 = diag(p)
  Sig1[1,2] = s
  Sig1[2,1] = s
  Sig1[3,4] = s
  Sig1[4,3] = s
  X1 = mvrnorm(n = kappa0[1], mu = rep(0,p), Sigma = Sig1)
  X2 = mvrnorm(n = kappa0[2]- kappa0[1], mu = c(rep(1,4), rep(0,4)), Sigma = Sig1)
  X3 = mvrnorm(n = N - kappa0[2], mu = c(rep(1,4), rep(0,4)), Sigma = diag(p))
  X = rbind(X1, X2, X3)
  x1 = X[,1]
  x2 = X[,2]
  x3 = X[,3]
  x4 = X[,4]
  x5 = X[,5]
  x6 = X[,6]
  x7 = X[,7]
  x8 = X[,8]
  x1_sq = X[,1]^2
  x2_sq = X[,2]^2
  x3_sq = X[,3]^2
  x4_sq = X[,4]^2
  x5_sq = X[,5]^2
  x6_sq = X[,6]^2
  x7_sq = X[,7]^2
  x8_sq = X[,8]^2
  x1x2 = X[,1]*X[,2]
  x1x3 = X[,1]*X[,3]
  x1x4 = X[,1]*X[,4]
  x1x5 = X[,1]*X[,5]
  x1x6 = X[,1]*X[,6]
  x1x7 = X[,1]*X[,7]
  x1x8 = X[,1]*X[,8]
  x2x3 = X[,2]*X[,3]
  x2x4 = X[,2]*X[,4]
  x2x5 = X[,2]*X[,5]
  x2x6 = X[,2]*X[,6]
  x2x7 = X[,2]*X[,7]
  x2x8 = X[,2]*X[,8]
  x3x4 = X[,3]*X[,4]
  x3x5 = X[,3]*X[,5]
  x3x6 = X[,3]*X[,6]
  x3x7 = X[,3]*X[,7]
  x3x8 = X[,3]*X[,8]
  x4x5 = X[,4]*X[,5]
  x4x6 = X[,4]*X[,6]
  x4x7 = X[,4]*X[,7]
  x4x8 = X[,4]*X[,8]
  x5x6 = X[,5]*X[,6]
  x5x7 = X[,5]*X[,7]
  x5x8 = X[,5]*X[,8]
  x6x7 = X[,6]*X[,7]
  x6x8 = X[,6]*X[,8]
  x7x8 = X[,7]*X[,8]
  # 2-degree polynomial feature embedding of X = (x1, x2, x3, x4)
  psi = cbind(x1, x2, x3, x4, x5, x6, x7, x8,
              x1_sq, x2_sq, x3_sq, x4_sq, x5_sq, x6_sq, x7_sq, x8_sq,
              x1x2, x1x3, x1x4, x1x5, x1x6, x1x7,x1x8,
              x2x3, x2x4, x2x5, x2x6, x2x7,x2x8,
              x3x4, x3x5, x3x6, x3x7, x3x8,
              x4x5, x4x6, x4x7, x4x8, 
              x5x6, x5x7, x5x8,
              x6x7, x6x8, 
              x7x8)
  data = list()
  data$X = X
  data$psi = psi
  return(data)
}

################## Generate independent series ###########################
kappa0 = c(100, 500)
N = 600
p = 8
psi_parallel <- function(i){
  sim_mvn_cicm(N = N, kappa0 = kappa0)
}
num_series = 100 # number of independent time series
data = mclapply(X = 1:num_series, FUN = psi_parallel, mc.cores = 10)
