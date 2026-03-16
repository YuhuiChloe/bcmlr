###################### A function to simulate psi ###########################
# A function that takes no inputs,
# simulates multivariate Gaussian X with TWO changes in covariances
# and returns the 2-degree polynomial embedding X
sim_mvn_cim <- function(N, kappa0, p){
  N = N
  kappa0 = kappa0
  p = p
  Sig1 = diag(p)
  X1 = mvrnorm(n = kappa0[1], mu = rep(0,p), Sigma = Sig1)
  X2 = mvrnorm(n = kappa0[2]- kappa0[1], mu = c(rep(2,2), rep(0,p-4), rep(-2,2)), Sigma = Sig1)
  X3 = mvrnorm(n = N - kappa0[2], mu = rep(0,p), Sigma = Sig1)
  X = rbind(X1, X2, X3)
  data = list()
  data$X = X # use a list to be consistent with CIC and CIMC functions
  return(data)
}

################## Generate independent series ###########################
N = 600 
kappa0 = c(100, 500)
p = 40
psi_parallel <- function(i){
  sim_mvn_cim(N = N, kappa0 = kappa0, p = p)
}
num_series = 100 # number of independent time series
data = mclapply(X = 1:num_series, FUN = psi_parallel, mc.cores = 10)
