
# Set the path
setwd("/Users/yuhuiwang/Desktop/ChangePoint/HighDim_Gaussian/Mispf-2trueCP")

#################### Data with 2 CPs ######################
# Read in data
X = readRDS("simdata_CPs155_255.rds")


########################## Call the bcmlr function #############################

# Fit 2 CP
num_iter = 10000
num_warmup = 2000
num_tune = 5000
num_CP = 2
init_num_temper = 40
out = bcmlr(data = X, num_iter = num_iter, num_warmup = num_warmup, num_temper = init_num_temper, 
            num_CP = num_CP, num_tune = num_tune, pc_cores = detectCores()/2)
saveRDS(out, file = "fit2CP.rds")


# Fit 3 CP
num_CP = 3
num_tune = 3000
out1 = bcmlr(data = X, num_iter = num_iter, num_warmup = num_warmup, num_temper = init_num_temper, 
            num_CP = num_CP, num_tune = num_tune, pc_cores = detectCores()/2)
saveRDS(out1, file = "fit3CP.rds")
