

#################### Data with 2 CPs ######################
# Set the path
setwd("/Users/yuhuiwang/Desktop/ChangePoint/HighDim_Gaussian/Mispf-2trueCP")
# Read in data
X = readRDS("simdata_CPs155_255.rds")


########################## Call the bcmlr function #############################

# Fit 2 CP
num_iter = 10000
num_warmup = 5000
num_tune = 5000
num_CP = 2
init_num_temper = 40
out = bcmlr(data = X, num_iter = num_iter, num_warmup = num_warmup, num_temper = init_num_temper, 
            num_CP = num_CP, num_tune = num_tune, pc_cores = detectCores()/2)
saveRDS(out, file = "fit2CP.rds")


# Fit 3 CP (Misspecified) 
num_CP = 3
out1 = bcmlr(data = X, num_iter = num_iter, num_warmup = num_warmup, num_temper = init_num_temper, 
            num_CP = num_CP, num_tune = num_tune, pc_cores = detectCores()/2)
saveRDS(out1, file = "fit3CP.rds")

#################### Data with 2 CPs ######################
# Set the path
setwd("/Users/yuhuiwang/Desktop/ChangePoint/HighDim_Gaussian/Mispf-1trueCP")
# Read in data
X = readRDS("simdata_CP_243.rds")

########################## Call the bcmlr function #############################

# Fit 1 CP
num_iter = 10000
num_warmup = 5000
num_tune = 5000 # not used when fitting 1 CP
num_CP = 1
init_num_temper = 40 # not used when fitting 1 CP
out = bcmlr(data = X, num_iter = num_iter, num_warmup = num_warmup, num_temper = init_num_temper, 
            num_CP = num_CP, num_tune = num_tune, pc_cores = detectCores()/2)
saveRDS(out, file = "fit1CP.rds")

# Fit 2 CP (Misspecified) 
num_CP = 2
out = bcmlr(data = X, num_iter = num_iter, num_warmup = num_warmup, num_temper = init_num_temper, 
            num_CP = num_CP, num_tune = num_tune, pc_cores = detectCores()/2)
saveRDS(out, file = "fit2CP.rds")

