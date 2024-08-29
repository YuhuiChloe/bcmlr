# data(wine)
# str(wine)
# 
# data(heartdisease)

N <- 400
kappa0 <- c(85, 170, 275)
p = 3
X1 = mvrnorm(n = kappa0[1], mu = runif(p, min = -2, max = 0), Sigma = diag(p))
X2 = mvrnorm(n = kappa0[2] - kappa0[1], mu = runif(p, min = 1, max = 2), Sigma = diag(p))
X3 = mvrnorm(n = kappa0[3] - kappa0[2], mu = runif(p, min = 3, max = 4), Sigma = diag(p))
X4 = mvrnorm(n = N - kappa0[3], mu = runif(p, min = 5, max = 6), Sigma = diag(p))

X = rbind(X1, X2, X3, X4)
X = scale(X, center = FALSE, scale = apply(X, 2, sd, na.rm = TRUE))
cov(X)
# Add an intercept:
X <- cbind(rep(1, N), X)
dim(X)




##### Use scatterplot3d #####
cluster <- factor(c(rep(1, nrow(X1)),
                    rep(2, nrow(X2)),
                    rep(3, nrow(X3)),
                    rep(4, nrow(X4))))
dim(X)
length(cluster)
library(scatterplot3d)
par(mfrow = c(1,1))
# Assuming X1, X2, X3, X4 are your four dimensions
scatterplot3d(X[,1], X[,2], X[,3], color=as.numeric(cluster),
              pch=19, xlab="Dimension 1", ylab="Dimension 2", zlab="Dimension 3",
              main="3D Scatter Plot of data X")
legend("topright", legend = paste("Cluster", 1:4), col = 1:4, pch = 19)

##### Use plotly ######
library(plotly)
cluster <- factor(c(rep(1, nrow(X1)),
                    rep(2, nrow(X2)),
                    rep(3, nrow(X3)),
                    rep(4, nrow(X4))))
# Create 3D scatter plot
fig <- plot_ly(x = X[,1], y = X[,2], z = X[,3], color = cluster, colors = c('#FF0000', '#00FF00', '#0000FF', '#FFFF00')) %>%
  add_markers() %>%
  layout(scene = list(xaxis = list(title = 'Dimension 1'),
                      yaxis = list(title = 'Dimension 2'),
                      zaxis = list(title = 'Dimension 3')),
         title = "3D Scatter Plot of data X")
print(fig)


X = as.matrix(X) 
dim(X)


######### Run the bcmlr function ########
out = bcmlr(data = X)



###### Results #######
runningtime = end - start
paste("The algorithm took ", format(runningtime), "to run ", 
      num_iter, "iterations, including ", num_warmup, "burn-in iterations.")


# Compute and print the mean coefficients (under tempering power 1)
mean_beta = apply(X=out$Beta, MARGIN=c(2,3), FUN=mean)
row_names <- paste("Beta", 1:nrow(mean_beta), sep = " ")
col_names <- paste("Class", 1:ncol(mean_beta), sep = " ")
dimnames(mean_beta) <- list(row_names, col_names)

kable(mean_beta, caption = "Posterior mean coefficients")


mean_Xbeta = X%*%mean_beta # Fitted values (linear) based on the posterior mean coefficients

mean_P = apply(out$P, MARGIN = c(2,3), FUN = mean)   # Posterior mean probabilities of success

J = 4
L = J-1
num_iter = 10000
num_warmup = 5000
num_temper = 30
for (j in 1:J){
  par(mfrow = c(3,1))
  plot(mean_P[,j], 
       main = paste("Class", j, ", Posterior mean probabilities of success"), 
       ylab="probabilities of success")
  plot(mean_Xbeta[,j], ylab = "Fitted values ", xlab = "Time points", 
       main = paste("Class", j, ", Fitted values based on the posterior mean coefficients"))
  if (j <= L){
    hist(out$Kappa[1:(num_iter-num_warmup),num_temper,j], breaks=100, xlim=c(0,N), 
         ylim = c(0,num_iter-num_warmup),
         ylab="Number of samples",
         main = paste("Histogram of #", j, "change point samples"))
    for (l in 1:length(kappa0)){
      abline(v=kappa0[l], col="red")
    }
  }
}

if (J == 4){
  Xbeta_plot <- plot_ly(x = mean_Xbeta[,1], y = mean_Xbeta[,2], z = mean_Xbeta[,3], 
                        color = cluster, 
                        colors = c('#FF0000', '#00FF00', '#0000FF', '#FFFF00')) %>%
    add_markers() %>%
    layout(scene = list(xaxis = list(title = 'Dimension 1'),
                        yaxis = list(title = 'Dimension 2'),
                        zaxis = list(title = 'Dimension 3')),
           title = "3D Scatter Plot of fitted values") 
  print(Xbeta_plot)
}


