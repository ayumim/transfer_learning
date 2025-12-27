library(mvtnorm)
library(invgamma)
source("gibbs.R")

p <- 10
n_t <- 300
sigma2 <- 1
# theta_s <- sample(c(runif(p, -6, -3),runif(p, 3, 6)), p)
theta_s <- ifelse(runif(p)<0.5, runif(p, -6, -3), runif(p, 3, 6))
x_t <- rmvnorm(n_t, rep(0, p), diag(1, p))
f_i <- c(x_t %*% theta_s)

kappa2 <- c(0.1, 10)

par(mfrow=c(1,3), oma=c(0,0,4,0))
results <- list()
for(i in 1:length(kappa2)) {
  theta_t <- rnorm(p, mean = theta_s, sd = sqrt(kappa2[i]))
  y <- rnorm(n_t, mean=x_t %*%theta_t, sd=sqrt(sigma2))
  results[[i]] <- gibbs(y, x_t, f_i)
  # results[[i]] <- burnin(results[[i]])
  plot(results[[i]]$delta, type="l")
  plot(results[[i]]$gamma, type="l")
  plot(results[[i]]$sigma, type="l")
  hist(results[[i]]$delta, xlab=expression(delta), main="")
  hist(results[[i]]$gamma, xlab=expression(gamma^2), main="")
  hist(results[[i]]$sigma, xlab=expression(sigma^2), main="")
  mtext(bquote(paste(kappa^2, "=", .(kappa2[i]))), line = 0, outer = TRUE)
}





