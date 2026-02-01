library(BayesLogit) # for rpg()
library(LearnBayes) # for computing inverse gamma: rigamma(), alternative to 1/rgamma()
library(MASS) # for mvrnorm()
library(parallel) # for parallel computing using mclapply() ---- Use this when you want to run bcmlr on independent data sets in parallel
library(pbmcapply) # to track and visualize the progress of parallel version of vectorized R functions
library(pdfCluster) # for computing Adjusted Rand Index
library(pracma) # for computing the Hausdorff distance
library(ecp) # The ecp packaeg includes E.divisive
library(reticulate) # for installing/loading Python packages
reticulate::py_require("ruptures") # Load ruptures for the current session 
rpt <- import("ruptures") # Import Python modules --- The ruptures package includes kcp and MultiRank. 
np <- import("numpy")
library(glmnet)
library(pROC) # for ROC curves and AUC values
options(pROCProgress = list(name = "none")) # Disable pROC progress bar globally

# In the bclr and bcmlr Gibbs samplers, we only need rpg() function to simulate the Polya Gamma variables. Rpg( ) is in BayesLogit. 
# We need to load PGdensity if we use parallel tempering. In the Metropolis-Hastings steps of parallel tempering, we need to compute the posterior, and that's where we need to use functions from PGdensity. 
library(PGdensity) # for PG density --- Users need to install this package following the steps in PGdensity/README.md. 
library(splines) # for monotone interpolation used in bcmlr-parallel tempering
