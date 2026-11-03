# Aula 8.2 - MCMC - Metropolis-Hastings e Gibbs - Exercicios

# 1 ----

## 1.1 ----

dens_gama = function(x, alfa, beta) {
  dgamma(x, shape = alfa, rate = beta)
}

curve(dens_gama(x, 2.5, 1), xlim = c(0, 10))

MCMC_MH_gama = function(M, tau2, x0 = 1, alfa = 2.5, beta = 1) {
  x = NULL
  x[1] = x0
  cont = 0
  
  for (t in 2:M) {
    x.p = rnorm(1, x[t-1], sqrt(tau2))
    
    if (x.p <= 0) {
      x[t] = x[t-1]
    } else {
      razao.num = dens_gama(x.p, alfa, beta)
      razao.den = dens_gama(x[t-1], alfa, beta)
      razao = razao.num / razao.den
      alpha = min(1, razao)
      u = runif(1, 0, 1)
      
      if (u < alpha) {
        x[t] = x.p
        cont = cont + 1
      } else {
        x[t] = x[t-1]
      }
    }
  }
  list(
    x = x,
    cont = cont
  )
}

q1 <- MCMC_MH_gama(20000, 1)

ts.plot(q1$x[(1001):10000])
acf(q1$x[(1001):10000])

hist(q1$x[seq(1001, 10000, by = 50)], freq = F)
curve(dens_gama(x, 2.5, 1), add = T)

q1$cont/20000

library(coda)
cadeia1 <- mcmc(q1$x)
geweke.diag(cadeia1)
effectiveSize(cadeia1)/20000

## 1.2 ----

q1_0.1 <- MCMC_MH_gama(10000, 0.1)
q1_1 <- MCMC_MH_gama(10000, 1)
q1_10 <- MCMC_MH_gama(10000, 10)
q1_100 <- MCMC_MH_gama(10000, 100)

q1_0.1$cont/10000
q1_1$cont/10000
q1_10$cont/10000
q1_100$cont/10000

cadeia_q1_0.1 <- mcmc(q1_0.1$x)
cadeia_q1_1 <- mcmc(q1_1$x)
cadeia_q1_10 <- mcmc(q1_10$x)
cadeia_q1_100 <- mcmc(q1_100$x)

geweke.diag(cadeia_q1_0.1)
geweke.diag(cadeia_q1_1)
geweke.diag(cadeia_q1_10)
geweke.diag(cadeia_q1_100)

effectiveSize(cadeia_q1_0.1)/10000
effectiveSize(cadeia_q1_1)/10000
effectiveSize(cadeia_q1_10)/10000
effectiveSize(cadeia_q1_100)/10000

# 2 ----

dens_mistura = function(x) {
  0.2 * dnorm(x, -7, sqrt(2)) +
  0.5 * dnorm(x,  0, sqrt(5)) +
  0.3 * dnorm(x,  4, 1)
}

curve(dens_mistura(x), xlim = c(-10,7))

MCMC_MH_mistura = function(M, sigma, x0 = 0) {
  x = NULL
  x[1] = x0
  cont = 0
  
  for (t in 2:M) {
    x.p = rnorm(1, x[t-1], sigma)
    
    razao = dens_mistura(x.p) / dens_mistura(x[t-1])
    alpha = min(1, razao)
    u = runif(1)
    
    if (u < alpha) {
      x[t] = x.p
      cont = cont + 1
    } else {
      x[t] = x[t-1]
    }
  }
  list(
    x = x,
    cont = cont
  )
}

q2_0.5 <- MCMC_MH_mistura(100000, 0.5)
q2_2 <- MCMC_MH_mistura(100000, 2)
q2_10 <- MCMC_MH_mistura(100000, 10)

q2_0.5$cont/100000
q2_2$cont/100000
q2_10$cont/100000

ts.plot(q2_0.5$x[10001:100000])
ts.plot(q2_2$x[10001:100000])
ts.plot(q2_10$x[10001:100000])

acf(q2_0.5$x[10001:100000])
acf(q2_2$x[10001:100000])
acf(q2_10$x[10001:100000])

hist(q2_0.5$x[10001:100000], freq = F)
hist(q2_2$x[10001:100000], freq = F)
hist(q2_10$x[10001:100000], freq = F)

# 3 ----
gibbs_exp_soma <- function(M, lambda, c, x0){
  n = length(lambda)
  x = matrix(NA, nrow = M, ncol = n)

  x[1,] = x0

  for(t in 2:M){
    x.p <- x[t-1, ]
    
    for(i in 1:n){
      soma = sum(x.p[-i])
      y = rexp(1, lambda[i])

      if(soma > c){
        x.p[i] = y
      } else {
        x.p[i] = (c - soma) + y
      }
    }

    x[t,] = x.p
  }
  x
}

q3 = gibbs_exp_soma(10000, c(1,2,3), 3, c(4,4,4))

ts.plot(q3[1001:10000, 1])
ts.plot(q3[1001:10000, 2])
ts.plot(q3[1001:10000, 3])

acf(q3p[, 1])
acf(q3p[, 2])
acf(q3p[, 3])

hist(q3[1001:10000, 1], freq = F)
hist(q3[1001:10000, 2], freq = F)
hist(q3[1001:10000, 3], freq = F)

# 4 ----

gibbs_normal <- function(M,y,mu0,tau0,a,b){
  n = length(y)
  ybar = mean(y)

  mu = NULL
  tau = NULL
  mu[1] = ybar
  tau[1] = 1/var(y)

  for(t in 2:M){

    media = (tau0*mu0 + n*tau[t-1]*ybar)/(tau0 + n*tau[t-1])
    variancia = 1/(tau0 + n*tau)
    mu[t] = rnorm(1, media, variancia)

    alfa = a + n/2
    beta = b + (1/2) * sum((y - mu[t])^2)
    tau[t] = rgamma(1, alfa, beta)
  }
  list(
    mu = mu,
    tau = tau
  )
}

n <- 100
mu_real <- 5
tau_real <- 1/4
y <- rnorm(n, mu_real, sqrt(1/tau_real))

q4 <- gibbs_normal(10000, y, 0, 0.01, 0.01, 0.01)

ts.plot(q4$mu[1001:10000])
ts.plot(q4$tau[1001:10000])

acf(q4$mu[1001:10000])
acf(q4$tau[1001:10000])

hist(q4$mu[1001:10000], freq = F)
abline(v = mu_real, col = "red")

hist(q4$tau[1001:10000], freq = F)
abline(v = tau_real, col = "red")
