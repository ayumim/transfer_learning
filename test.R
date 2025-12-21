library(mvtnorm)
library(invgamma)

p <- 10
n_t <- 100
sigma2 <- 1
# theta_s <- sample(c(runif(p, -6, -3),runif(p, 3, 6)), p)
theta_s <- ifelse(runif(p)<0.5, runif(p, -6, -3), runif(p, 3, 6))

kappa2 <- 0.1
theta_t <- rnorm(p, mean = theta_s, sd = sqrt(kappa2))

x_t <- rmvnorm(n_t, rep(0, p), diag(1, p))
y <- rnorm(n_t, mean=x_t %*%theta_t, sd=sqrt(sigma2))

f_i <- c(x_t %*% theta_s)

gibbs <- function(y, x, f_i, N=10000, tau2=1, a_gamma=2, b_gamma=2, 
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
      alpha[i, j] <- rnorm(1, mean=m * v, sd=sqrt(v))
    }
    
    # delta
    v_delta <- 1 / (n / gamma_prev + 1 / tau2)
    m_delta <- (sum(alpha[i,])/gamma_prev) * v_delta 
    delta[i] <- rnorm(1, mean=m_delta, sd=sqrt(v_delta))
    
    # gamma
    shape_g <- a_gamma + n / 2
    rate_g <- sum((alpha[i,] - delta[i])^2) / 2 + b_gamma
    gamma[i] <- rinvgamma(1, shape=shape_g, rate=rate_g)
    
    # sigma 
    shape_s <- a_sigma + n / 2
    rate_s <- sum((y - (alpha[i,] + f_i))^2) / 2 + b_sigma
    sigma[i] <- rinvgamma(1, shape = shape_s, rate=rate_s)
  }
  
  return(list(alpha=alpha, delta=delta, gamma=gamma, sigma=sigma))
}
