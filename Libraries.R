library(BayesLogit) # for rpg()
library(MASS) # for mvrnorm()
library(PGdensity) # for PG density (users need to pre-download this package)
library(parallel) # for parallel computing using mclapply()
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
library(splines) # for monotone interpolation used in bcmlr-parallel tempering
