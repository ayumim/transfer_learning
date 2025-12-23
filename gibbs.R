# gibbs sampling functions 
gibbs <- function(y, x, f_i, N=10000, tau2=1, a_gamma=2, b_gamma=2, 
                  a_sigma=2, b_sigma=2) {
  n <- length(y)
  p <- ncol(x)
  
  eta <- matrix(NA, nrow=N, ncol=n)
  gamma <- numeric(N)
  delta <- numeric(N)
  sigma <- numeric(N)
  # initial value
  delta[1] <- 0
  gamma[1] <- 1
  sigma[1] <- 1
  eta[1,] <- rep(0, n)
  
  # start gibbs sampling
  for(i in 2:N) {
    eta_prev <- eta[i-1,]
    gamma_prev <- gamma[i-1]
    delta_prev <- delta[i-1]
    sigma_prev <- sigma[i-1]
    # eta
    for(j in 1:n) {
      # v <- 1 / (1/sigma_prev + 1/gamma_prev)
      # m <- (y[j] - f_i[j]) / sigma_prev + delta_prev / gamma_prev
      eta[i, j] <- eta_update(y[j], f_i[j], sigma_prev, gamma_prev, delta_prev)
    } # end of j
    
    # delta
    delta[i] <- delta_update(eta[i,], gamma_prev, tau2, n)
    
    # v_delta <- 1 / (n / gamma_prev + 1 / tau2)
    # m_delta <- (sum(eta[i,])/gamma_prev) * v_delta 
    # delta[i] <- rnorm(1, mean=m_delta, sd=sqrt(v_delta))
    
    # gamma
    gamma[i] <- gamma_update(eta[i,], delta[i], a_gamma, b_gamma, n)
    # shape_g <- a_gamma + n / 2
    # rate_g <- sum((eta[i,] - delta[i])^2) / 2 + b_gamma
    # gamma[i] <- rinvgamma(1, shape=shape_g, rate=rate_g)
    
    # sigma 
    sigma[i] <- sigma_update(y, eta[i,], f_i, a_sigma, b_sigma, n)
    # shape_s <- a_sigma + n / 2
    # rate_s <- sum((y - (eta[i,] + f_i))^2) / 2 + b_sigma
    # sigma[i] <- rinvgamma(1, shape = shape_s, rate=rate_s)
  } # end of i
  
  return(list(eta=eta, delta=delta, gamma=gamma, sigma=sigma))
}

# eta update 
eta_update <- function(y, f_i, sigma, gamma, delta) {
  v <- 1 / (1/sigma + 1/gamma)
  m <- (y - f_i) / sigma + delta / gamma
  return(rnorm(1, mean=m * v, sd=sqrt(v)))
}

# delta update
delta_update <- function(eta, gamma, tau2, n) {
  v <- 1 / (n / gamma + 1 / tau2)
  m <- (sum(eta) / gamma) * v
  return(rnorm(1, mean=m, sd=sqrt(v)))
}

# gamma update 
gamma_update <- function(eta, delta, a, b, n) {
  shape <- a + n / 2
  rate <- sum((eta - delta)^2) / 2 + b
  return(rinvgamma(1, shape=shape, rate=rate))
}

# sigma update 
sigma_update <- function(y, eta, f_i, a, b, n) {
  shape <- a + n / 2
  rate <- sum((y - eta - f_i)^2) / 2 + b
  return(rinvgamma(1, shape = shape, rate=rate))
}