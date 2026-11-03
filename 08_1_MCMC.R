# Aula 8.1 - MCMC - Metropolis-Hastings e Gibbs 

library(coda)
library(mvtnorm)

# Metropolis-Hastings ----

## 1.1 ----
# Inversa-Gama(a, b) com proposta Uniforme(0, L)

dens_inv_gama <- function(x, a, b) {
  ifelse(x > 0, (b^a / gamma(a)) * x^(-a - 1) * exp(-b / x), 0)
}

curve(dens_inv_gama(x, 2, 2), 0, 10)

MCMC_MH_InvGama_Uniforme <- function(M, x0, a, b, L = 100) {
  x <- numeric(M)
  x[1] <- x0
  cont <- 0
  
  for (t in 2:M) {
    x.p <- runif(1, 0, L)
    razao <- dens_inv_gama(x.p, a, b) / dens_inv_gama(x[t-1], a, b)
    alpha <- min(1, razao)
    u <- runif(1)
    
    if (u < alpha) {
      x[t] <- x.p
      cont <- cont + 1
    } else {
      x[t] <- x[t-1]
    }
  }
  list(
    x = x,
    taxa_aceitacao = cont / M
  )
}

q1_unif <- MCMC_MH_InvGama_Uniforme(10000, 0.1, 2, 3)
q1_unif$taxa_aceitacao

ts.plot(q1_unif$x[1001:10000])
acf(q1_unif$x[1001:10000])
hist(q1_unif$x[seq(1001, 10000, by = 40)], freq = FALSE, breaks = 50)
curve(dens_inv_gama(x, 2, 3), col = "red", add = TRUE)

## 1.2 ----
# Inversa-Gama(a, b) com proposta Gama(shape_prop, rate_prop)

MCMC_MH_InvGama_Gama <- function(M, x0, a, b, shape_prop, rate_prop) {
  x <- numeric(M)
  x[1] <- x0
  cont <- 0
  
  for (t in 2:M) {
    x.p <- rgamma(1, shape = shape_prop, rate = rate_prop)
    
    # Razao de Hastings completa: (f(xp)/f(xc)) * (q(xc|xp)/q(xp|xc))
    num <- dens_inv_gama(x.p, a, b) * dgamma(x[t-1], shape_prop, rate_prop)
    den <- dens_inv_gama(x[t-1], a, b) * dgamma(x.p, shape_prop, rate_prop)
    razao <- num / den
    alpha <- min(1, razao)
    u <- runif(1)
    
    if (u < alpha) {
      x[t] <- x.p
      cont <- cont + 1
    } else {
      x[t] <- x[t-1]
    }
  }
  list(
    x = x,
    taxa_aceitacao = cont / M
  )
}

q1_gama <- MCMC_MH_InvGama_Gama(10000, 0.1, 2, 3, 2/3, 1/3)
q1_gama$taxa_aceitacao

ts.plot(q1_gama$x[1001:10000])
acf(q1_gama$x[1001:10000])
hist(q1_gama$x[seq(1001, 10000, by = 20)], freq = FALSE, breaks = 50)
curve(dens_inv_gama(x, 2, 3), col = "red", add = TRUE)

cadeia1 <- mcmc(q1_gama$x, start = 1001, thin = 30)
summary(cadeia1)
geweke.diag(cadeia1)   # deve estar entre -1.96 e 1.96
effectiveSize(cadeia1)/10000

## 2 ----
# t de Student com proposta Normal(mu, 1)

dt3p <- function(x, mu, sigma2, nu) {
  constante <- gamma((nu + 1) / 2) / (gamma(nu / 2) * sqrt(pi * nu * sigma2))
  nucleo <- (1 + (1 / nu) * ((x - mu) / sqrt(sigma2))^2)^(-(nu + 1) / 2)
  return(constante * nucleo)
}

mu <- 5
sigma2 <- 1
nu <- 10
curve(dt3p(x, mu, sigma2, nu), 0, 10)

MCMC_MH_t3p_Normal <- function(M, mu, sigma2, nu, x0 = 0.1) {
  x <- numeric(M)
  x[1] <- x0
  cont <- 0
  
  for (t in 2:M) {
    x.p <- rnorm(1, mean = mu, sd = 1)
    
    num <- dt3p(x.p, mu, sigma2, nu) * dnorm(x[t-1], mean = mu, sd = 1)
    den <- dt3p(x[t-1], mu, sigma2, nu) * dnorm(x.p, mean = mu, sd = 1)
    razao <- num / den
    alpha <- min(1, razao)
    u <- runif(1)
    
    if (u < alpha) {
      x[t] <- x.p
      cont <- cont + 1
    } else {
      x[t] <- x[t-1]
    }
  }
  list(
    x = x,
    taxa_aceitacao = cont / M
  )
}

q2_normal <- MCMC_MH_t3p_Normal(M = 10000, mu = mu, sigma2 = sigma2, nu = nu)
q2_normal$taxa_aceitacao

ts.plot(q2_normal$x[1001:10000])
hist(q2_normal$x[1001:10000], freq = FALSE, breaks = 50)
curve(dt3p(x, mu, sigma2, nu), col = "red", add = TRUE)

cadeia2 <- mcmc(q2_normal$x, start = 1001, thin = 30)
summary(cadeia2)
geweke.diag(cadeia2)
effectiveSize(cadeia2)/10000

## 3 ----
# t de Student com proposta N(x[t-1], tau2) PASSEIO ALEATORIO

MCMC_MH_t3p_RW <- function(M, mu, sigma2, nu, tau2, x0 = 0.1) {
  x <- numeric(M)
  x[1] <- x0
  cont <- 0
  
  for (t in 2:M) {
    x.p <- rnorm(1, mean = x[t-1], sd = sqrt(tau2))
    
    razao <- dt3p(x.p, mu, sigma2, nu) / dt3p(x[t-1], mu, sigma2, nu)
    alpha <- min(1, razao)
    u <- runif(1)
    
    if (u < alpha) {
      x[t] <- x.p
      cont <- cont + 1
    } else {
      x[t] <- x[t-1]
    }
  }
  list(
    x = x,
    taxa_aceitacao = cont / M
  )
}

q2_rw_01 <- MCMC_MH_t3p_RW(M = 10000, mu = mu, sigma2 = sigma2, nu = nu, tau2 = 0.1)
q2_rw_1  <- MCMC_MH_t3p_RW(M = 10000, mu = mu, sigma2 = sigma2, nu = nu, tau2 = 1)
q2_rw_10 <- MCMC_MH_t3p_RW(M = 10000, mu = mu, sigma2 = sigma2, nu = nu, tau2 = 10)

q2_rw_01$taxa_aceitacao
q2_rw_1$taxa_aceitacao
q2_rw_10$taxa_aceitacao

ts.plot(q2_rw_01$x[1001:10000])
ts.plot(q2_rw_1$x[1001:10000])
ts.plot(q2_rw_10$x[1001:10000])

## 4 ----
# mistura de duas normais bivariadas com pesos pi_mix e 1 - pi_mix
# PASSEIO ALEATORIO

mu1 <- c(0, 6)
Sigma1 <- matrix(c(4, 0.5, 0.5, 2), nrow = 2)
mu2 <- c(-5, 4)
Sigma2 <- matrix(c(1, -0.2, -0.2, 1), nrow = 2)
pi_mix <- 0.4

dens_mistura_biv <- function(vec, mu1, Sigma1, mu2, Sigma2, pi_mix) {
  pi_mix * dmvnorm(vec, mean = mu1, sigma = Sigma1) +
    (1 - pi_mix) * dmvnorm(vec, mean = mu2, sigma = Sigma2)
}

MCMC_MH_MisturaBiv <- function(M, mu1, Sigma1, mu2, Sigma2, pi_mix, sigma_prop = 1, x0 = c(0, 0)) {
  x_vec <- matrix(NA, nrow = 2, ncol = M)
  x_vec[, 1] <- x0
  cont <- 0
  
  for (t in 2:M) {
    x.p <- rmvnorm(1, mean = x_vec[, t-1], sigma = diag(sigma_prop^2, 2))
    
    razao <- dens_mistura_biv(x.p, mu1, Sigma1, mu2, Sigma2, pi_mix) /
             dens_mistura_biv(x_vec[, t-1], mu1, Sigma1, mu2, Sigma2, pi_mix)
    alpha <- min(1, razao)
    u <- runif(1)
    
    if (u < alpha) {
      x_vec[, t] <- x.p
      cont <- cont + 1
    } else {
      x_vec[, t] <- x_vec[, t-1]
    }
  }
  list(
    x_vec = x_vec,
    taxa_aceitacao = cont / M
  )
}


q4 <- MCMC_MH_MisturaBiv(M = 50000, mu1 = mu1, Sigma1 = Sigma1, mu2 = mu2, Sigma2 = Sigma2, pi_mix = pi_mix, sigma_prop = 1.5)
q4$taxa_aceitacao

ts.plot(q4$x_vec[1, 10001:50000])
ts.plot(q4$x_vec[2, 10001:50000])

cadeia4 <- mcmc(t(q4$x_vec), start = 10001, thin = 50)
summary(cadeia4)
geweke.diag(cadeia4)
effectiveSize(cadeia4)/10000

# Gibbs ----

## 5.1 ----
# Normal bivariada padrao 
# X | Y = y ~ N(rho * y, 1 - rho^2)
# Y | X = x ~ N(rho * x, 1 - rho^2).

Gibbs_NormalBiv <- function(M, rho, x0 = 10, y0 = -5) {
  x <- numeric(M)
  y <- numeric(M)
  x[1] <- x0
  y[1] <- y0
  
  sd_cond <- sqrt(1 - rho^2)
  
  for (t in 2:M) {
    # Atualiza X dado Y[t-1]
    x[t] <- rnorm(1, mean = rho * y[t-1], sd = sd_cond)
    # Atualiza Y dado X[t] (recem-atualizado, caracteristica do Gibbs)
    y[t] <- rnorm(1, mean = rho * x[t], sd = sd_cond)
  }
  list(
    x = x,
    y = y
  )
}

rho <- -0.8
q5 <- Gibbs_NormalBiv(M = 10000, rho = rho)

ts.plot(q5$x[1001:10000])
ts.plot(q5$y[1001:10000])

cadeia5 <- mcmc(cbind(q5$x, q5$y), start = 1001, thin = 20)
summary(cadeia5)
geweke.diag(cadeia5)
effectiveSize(cadeia5)

## 6 ----
# Gibbs com passos de Metropolis-Hastings
# Mistura de normais bivariadas

Gibbs_MH_MisturaBiv <- function(M, mu1, Sigma1, mu2, Sigma2, pi_mix, tau2 = 1, x0 = c(0, 0)) {
  x_vec <- matrix(NA, nrow = 2, ncol = M)
  x_vec[, 1] <- x0
  cont1 <- 0
  cont2 <- 0
  
  for (t in 2:M) {
    # Atualiza X1 dado X2[t-1]
    x1.p <- rnorm(1, mean = x_vec[1, t-1], sd = sqrt(tau2))
    
    razao <- dens_mistura_biv(c(x1.p, x_vec[2, t-1]), mu1, Sigma1, mu2, Sigma2, pi_mix) /
             dens_mistura_biv(c(x_vec[1, t-1], x_vec[2, t-1]), mu1, Sigma1, mu2, Sigma2, pi_mix)
    alpha <- min(1, razao)
    u <- runif(1)
    
    if (u < alpha) {
      x_vec[1, t] <- x1.p
      cont1 <- cont1 + 1
    } else {
      x_vec[1, t] <- x_vec[1, t-1]
    }
    
    # Atualiza X2 dado X1[t] (recem-atualizado)
    x2.p <- rnorm(1, mean = x_vec[2, t-1], sd = sqrt(tau2))
    
    razao <- dens_mistura_biv(c(x_vec[1, t], x2.p), mu1, Sigma1, mu2, Sigma2, pi_mix) /
             dens_mistura_biv(c(x_vec[1, t], x_vec[2, t-1]), mu1, Sigma1, mu2, Sigma2, pi_mix)
    alpha <- min(1, razao)
    u <- runif(1)
    
    if (u < alpha) {
      x_vec[2, t] <- x2.p
      cont2 <- cont2 + 1
    } else {
      x_vec[2, t] <- x_vec[2, t-1]
    }
  }
  list(
    x_vec = x_vec,
    taxa_aceitacao_x1 = cont1 / M,
    taxa_aceitacao_x2 = cont2 / M
  )
}

q6 <- Gibbs_MH_MisturaBiv(M = 50000, mu1 = mu1, Sigma1 = Sigma1, mu2 = mu2, Sigma2 = Sigma2, pi_mix = pi_mix, tau2 = 1)
q6$taxa_aceitacao_x1
q6$taxa_aceitacao_x2

ts.plot(q6$x_vec[5000, 5000:50000])
ts.plot(q6$x_vec[5000, 5000:50000])

hist(q6$x_vec[1, 5000:50000], freq = FALSE, breaks = 50)
hist(q6$x_vec[2, 5000:50000], freq = FALSE, breaks = 50)

cadeia6 <- mcmc(t(q6$x_vec), start = 5000, thin = 50)
summary(cadeia6)
geweke.diag(cadeia6)
effectiveSize(cadeia6)