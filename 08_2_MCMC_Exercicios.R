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

q1 <- MCMC_MH_gama(10000, 1)
q1$x

ts.plot(q1$x[(1001):10000])
acf(q1$x[(1001):10000])

hist(x[seq(1001, 10000, by = 50)], freq = F)
curve(dens_gama(x, 2.5, 1), add = T)

q1$cont/10000

## 1.2 ----

q1_0.1 <- MCMC_MH_gama(10000, 0.1)
q1_1 <- MCMC_MH_gama(10000, 1)
q1_10 <- MCMC_MH_gama(10000, 10)
q1_100 <- MCMC_MH_gama(10000, 100)

q1_0.1$cont/10000
q1_1$cont/10000
q1_10$cont/10000
q1_100$cont/10000

# 2 ----

dens_mistura = function(x) {
  0.2 * dnorm(x, -7, sqrt(2)) +
  0.5 * dnorm(x,  0, sqrt(5)) +
  0.3 * dnorm(x,  4, 1)
}

x = seq(-15, 10, by = 0.1)
y = dens_mistura(x)

plot(x, y, type = "l")

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

par(mfrow = c(3, 1))
ts.plot(q2_0.5$x[10001:100000])
ts.plot(q2_2$x[10001:100000])
ts.plot(q2_10$x[10001:100000])
par(mfrow = c(1, 1))

par(mfrow = c(3, 1))
acf(q2_0.5$x[10001:100000])
acf(q2_2$x[10001:100000])
acf(q2_10$x[10001:100000])
par(mfrow = c(1, 1))

par(mfrow = c(3,1))
for (i in 1:3) {
  q2 = list(q2_0.5, q2_2, q2_10)[[i]]
  hist(q2$x[(1001):100000], freq = F)
  curve(dens_mistura(x), add = T, col = "red")
}
par(mfrow = c(1, 1))

# 3 ----

gibbs_exp_soma = function(M, lambda, c, x0) {
  n = length(lambda)
  X = matrix(NA, nrow = M, ncol = n)
  
  X[1, ] = x0
  
  for (t in 2:M) {
    x.atual = X[t-1, ]
    
    for(i in 1:n){
      i = sample(1:n, 1)
      soma_outros = sum(x.atual[-i])

      #############################################################
      
      X[t, ] = x.atual
    }
  }
  X
}

# 4 ----
# Y_i ~ N(mu, 1/tau)
#
# mu ~ N(mu0, 1/tau0)
# tau ~ Gama(a, b)
#
# mu | tau, y ~ N( (tau0*mu0 + n*tau*ybar) / (tau0 + n*tau),  1/(tau0 + n*tau) )
# tau | mu, y ~ Gama( a + n/2,  b + 0.5 * sum((y - mu)^2) )   (taxa b)

## 4.1 ----
# Dados simulados (valores "verdadeiros": mu = 5, sigma = 2 -> tau = 0.25)
 
n = 100
mu_real = 5
tau_real = 1/4
y = rnorm(n, mu_real, sqrt(1/tau_real))
 
ybar = mean(y)
 
# Hiperparametros das prioris (pouco informativas)
mu0 = 0
tau0 = 0.01
a = 0.01
b = 0.01
 
## 4.2 ----
# Gibbs
 
M = 10000
burn = 1000
 
mu = numeric(M)
tau = numeric(M)
 
# Valores iniciais
mu[1] = ybar
tau[1] = 1 / var(y)
 
for (t in 2:M) {
  # mu | tau, y
  prec = tau0 + n * tau[t-1]
  media = (tau0 * mu0 + n * tau[t-1] * ybar) / prec
  mu[t] = rnorm(1, media, sqrt(1 / prec))
  
  # tau | mu, y
  tau[t] = rgamma(1, shape = a + n/2, rate = b + 0.5 * sum((y - mu[t])^2))
}
 
## 4.3 ----
# Diagnosticos e resultados
 
par(mfrow = c(2, 1))
ts.plot(mu[(burn+1):M],  xlab = "iteracoes", ylab = expression(mu))
ts.plot(tau[(burn+1):M], xlab = "iteracoes", ylab = expression(tau))
par(mfrow = c(1, 1))
 
par(mfrow = c(1, 2))
acf(mu[(burn+1):M],  main = expression(mu))
acf(tau[(burn+1):M], main = expression(tau))
par(mfrow = c(1, 1))
 
par(mfrow = c(1, 2))
hist(mu[(burn+1):M],  freq = F, main = expression(mu),  xlab = expression(mu))
abline(v = mu_real, col = "red", lwd = 2)
hist(tau[(burn+1):M], freq = F, main = expression(tau), xlab = expression(tau))
abline(v = tau_real, col = "red", lwd = 2)
par(mfrow = c(1, 1))
 
# Resumos a posteriori: media e intervalo de credibilidade de 95%
mu_post  = mu[(burn+1):M]
tau_post = tau[(burn+1):M]
 
rbind(mu  = c(media = mean(mu_post),  quantile(mu_post,  c(0.025, 0.975))),
      tau = c(media = mean(tau_post), quantile(tau_post, c(0.025, 0.975))))
 
# Em termos da variancia sigma^2 = 1/tau
sigma2_post = 1 / tau_post
c(media = mean(sigma2_post), quantile(sigma2_post, c(0.025, 0.975)))
 
# Distribuicao conjunta
plot(mu_post, tau_post, pch = 16, cex = 0.3, xlab = expression(mu),
     ylab = expression(tau))
 