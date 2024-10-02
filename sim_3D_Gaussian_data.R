####################### a 3D data set with 1 CP #########################
kappa0 = 243
N = 400
p = 3
X1 = mvrnorm(n = kappa0, mu = runif(p, min = -1, max = 0.1), Sigma = diag(p))
X2 = mvrnorm(n = N - kappa0, mu = runif(p, min = 0.5, max = 2), Sigma = diag(p))
X = rbind(X1, X2)
X = scale(X, center = FALSE, scale = apply(X, 2, sd, na.rm = TRUE))
cov(X)
# add an intercept
X <- cbind(rep(1, N), X) 
dim(X)

cluster <- factor(c(rep(1, nrow(X1)),
                    rep(2, nrow(X2))))
# Use scatterplot3d #
length(cluster)
library(scatterplot3d)
par(mfrow = c(1,1))
scatterplot3d(X[,2], X[,3], X[,4], color=as.numeric(cluster),
              pch=19, xlab="Dimension 1", ylab="Dimension 2", zlab="Dimension 3",
              main="3D Scatter Plot of data X")
legend("topright", legend = paste("Cluster", 1:(length(kappa0)+1)), col = 1:(length(kappa0)+1), pch = 19)

# Set the path and store the data 
setwd("/Users/yuhuiwang/Desktop/ChangePoint/HighDim_Gaussian/Mispf-1trueCP")
saveRDS(X, file = "simdata_CP_243.rds") 





######## a 3D data set with 2 CPs##########
N <- 400
kappa0 <- c(155, 255)
p = 3
X1 = mvrnorm(n = kappa0[1], mu = runif(p, min = -1, max = 0), Sigma = diag(p))
X2 = mvrnorm(n = kappa0[2] - kappa0[1], mu = runif(p, min = 0.5, max = 2), Sigma = diag(p))
X3 = mvrnorm(n = N - kappa0[2], mu = runif(p, min = 2, max = 3), Sigma = diag(p))

X = rbind(X1, X2, X3)
X = scale(X, center = FALSE, scale = apply(X, 2, sd, na.rm = TRUE))
cov(X)
# add an intercept
X <- cbind(rep(1, N), X) 
dim(X)

# 3D Data Visualization
cluster <- factor(c(rep(1, nrow(X1)),
                    rep(2, nrow(X2)),
                    rep(3, nrow(X3))))
library(scatterplot3d)
par(mfrow = c(1,1))
scatterplot3d(X[,2], X[,3], X[,4], color=as.numeric(cluster),
              pch=19, xlab="Dimension 1", ylab="Dimension 2", zlab="Dimension 3",
              main="3D Scatter Plot of data X")
legend("topright", legend = paste("Cluster", 1:(length(kappa0)+1)), col = 1:(length(kappa0)+1), pch = 19)


# Set the path and store the data 
setwd("/Users/yuhuiwang/Desktop/ChangePoint/HighDim_Gaussian/Mispf-2trueCP")
saveRDS(X, file = "simdata_CPs155_255.rds") 


