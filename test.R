library(mvtnorm)

p <- 10
n_t <- 10
sigma2 <- 1
# theta_s <- sample(c(runif(p, -6, -3),runif(p, 3, 6)), p)
theta_s <- ifelse(runif(p)<0.5, runif(p, -6, -3), runif(p, 3, 6))

kappa2 <- diag(10, p)
theta_t <- rnorm(p, mean = theta_s, sd = sqrt(kappa2))

x_t <- rmvnorm(n_t, rep(0, p), diag(1, p))



y <- rnorm(n_t, mean=x_t %*%theta_t, sd=sqrt(sigma2))

f_i <- c(x_t %*% theta_s)

gibbs <- function(y, x, N=10000, tau2=1, a_gamma=2, b_gamma=2, 
                  a_sigma=2, b_sigma=2) {
  n <- length(y)
  p <- ncol(x)
  
  alpha <- matrix(NA, nrow=N, ncol=n)
  gamma <- numeric(N)
  delta <- numeric(N)
  sigma <- numeric(N)
  # initial value
  delta[1] <- 0
  gamma[1] <- 1
  sigma[1] <- 1
  alpha[1,] <- rep(0, n)
  
  # gibbs 
  for(i in 2:N) {
    alpha_prev <- alpha[i-1,]
    gamma_prev <- gamma[i-1]
    delta_prev <- delta[i-1]
    sigma_prev <- sigma[i-1]
    # alpha
    for(j in 1:n) {
      v <- 1 / (1/sigma_prev + 1/gamma_prev)
      m <- (y[j] - f_i[j]) / sigma_prev + delta_prev / gamma_prev
      alpha[i-1, j] <- rnorm(1, mean=m * v, sd=sqrt(v))
    }
    
  }
}