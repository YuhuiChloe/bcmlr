# A generalized Bayesian approach to multiple changepoint analysis

The bcmlr function gives the options to 
- 1) detect a fixed number of changepoints when the number of changes is known OR
- 2) detect the number of changes when the number of changes is unknown. When the number of changes is unknown, one needs to use the bcmlr_model_select function that fits bcmlr to first detect the number of changes and generates posterior samples for posterior inference. 

# Main arguments in the bcmlr or bcmlr_model_select functions: 
- data: A data frame or matrix with numerical entries
- num_CP: The total number of changepoints the user anticipates (Notation in the paper: $k$)
- init: initialization method (The default initialization is to place the changepoints at evenly spaced positions of the data series.)
- prior: Choose "Gaussian" if you have low-dimensional data or "horseshoe" if you have high-dimensional data 
- alpha_f: default to 0.1. (1-alpha_f) is the confidence level o a frequentist confidence internval (Notation in the paper: $\alpha$)
- tempering: default to 1. It is only used when implementing non-reversible parallel tempering (Notation in the paper supplement: $t$)
- thinning: default at 1. It leads to how much data is held out (Notation in the paper:$\zeta)

# Outputs of the bcmlr_model_select function: 
- Kappa: posterior modes of the changepoint distribution (Notation in the paper: $\kappa$)
- Beta: posterior coefficients (from logistic regression or multinomial logistic regression) (Notation in the paper: $\beta$)
- P: posterior probablities (from logistic regression or multinomial logistic regression) (Notation in the paper: $p$)
- num_cp_dist: posterior distribution of the number of changepoints
- num_cp: the posterior mode of num_cp_dist = the length of Kappa (Notation in the paper: $L$ (truth), $\hat{k}$ (estimate))
- max_cps: the posterior modes based on the max number of changepoints the user anticipates (Notation in the paper: $k$)
- probs_num_cp: posterior probabilities of all possible numbers of changepoints (up to the length of max_cps). The number with the highest posterior probabilities is the number of changes bcmlr deem there to be. 

# Steps for single/multiple changepoint detection. 
- Download the PGDensity package
- Use the bcmlr function in bcmlr.R if the true number of changepoints is known. 
- Use the bcmlr_model_select function in bcmlr_model_select.R if the true number of changepoints is unknown. 


> [!IMPORTANT]
> This repository is still under construction.
