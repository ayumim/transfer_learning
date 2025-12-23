# gibbs sampling functions 
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
  
  # start gibbs sampling
  for(i in 2:N) {
    alpha_prev <- alpha[i-1,]
    gamma_prev <- gamma[i-1]
    delta_prev <- delta[i-1]
    sigma_prev <- sigma[i-1]
    # alpha
    for(j in 1:n) {
      # v <- 1 / (1/sigma_prev + 1/gamma_prev)
      # m <- (y[j] - f_i[j]) / sigma_prev + delta_prev / gamma_prev
      alpha[i, j] <- alpha_update(y[j], f_i[j], sigma_prev, gamma_prev, delta_prev)
    } # end of j
    
    # delta
    delta[i] <- delta_update(alpha[i,], gamma_prev, tau2, n)
    
    # v_delta <- 1 / (n / gamma_prev + 1 / tau2)
    # m_delta <- (sum(alpha[i,])/gamma_prev) * v_delta 
    # delta[i] <- rnorm(1, mean=m_delta, sd=sqrt(v_delta))
    
    # gamma
    gamma[i] <- gamma_update(alpha[i,], delta[i], a_gamma, b_gamma, n)
    # shape_g <- a_gamma + n / 2
    # rate_g <- sum((alpha[i,] - delta[i])^2) / 2 + b_gamma
    # gamma[i] <- rinvgamma(1, shape=shape_g, rate=rate_g)
    
    # sigma 
    sigma[i] <- sigma_update(y, alpha[i,], f_i, a_sigma, b_sigma, n)
    # shape_s <- a_sigma + n / 2
    # rate_s <- sum((y - (alpha[i,] + f_i))^2) / 2 + b_sigma
    # sigma[i] <- rinvgamma(1, shape = shape_s, rate=rate_s)
  } # end of i
  
  return(list(alpha=alpha, delta=delta, gamma=gamma, sigma=sigma))
}

# alpha update 
alpha_update <- function(y, f_i, sigma, gamma, delta) {
  v <- 1 / (1/sigma + 1/gamma)
  m <- (y - f_i) / sigma + delta / gamma
  return(rnorm(1, mean=m * v, sd=sqrt(v)))
}

# delta update
delta_update <- function(alpha, gamma, tau2, n) {
  v <- 1 / (n / gamma + 1 / tau2)
  m <- (sum(alpha) / gamma) * v
  return(rnorm(1, mean=m, sd=sqrt(v)))
}

# gamma update 
gamma_update <- function(alpha, delta, a, b, n) {
  shape <- a + n / 2
  rate <- sum((alpha - delta)^2) / 2 + b
  return(rinvgamma(1, shape=shape, rate=rate))
}

# sigma update 
sigma_update <- function(y, alpha, f_i, a, b, n) {
  shape <- a + n / 2
  rate <- sum((y - alpha - f_i)^2) / 2 + b
  return(rinvgamma(1, shape = shape, rate=rate))
}