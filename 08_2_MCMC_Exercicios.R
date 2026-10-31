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
    
    for (i in 1:n) {
      soma_outros = sum(x.atual[-i])
      y = rexp(1, rate = lambda[i])
      
      if (soma_outros > c) {
        x.atual[i] = y
      } else {
        x.atual[i] = (c - soma_outros) + y
      }
    }
    
    X[t, ] = x.atual
  }
  X
}

q3 = gibbs_exp_soma(
  10000,
  c(1,2,3),
  3,
  c(4,4,4)
)

par(mfrow = c(3, 1))
ts.plot(q3[1001:10000, 1], ylab = expression(x[1]))
ts.plot(q3[1001:10000, 2], ylab = expression(x[2]))
ts.plot(q3[1001:10000, 3], ylab = expression(x[3]))
par(mfrow = c(1, 1))

par(mfrow = c(3, 1))
acf(q3p[, 1])
acf(q3p[, 2])
acf(q3p[, 3])
par(mfrow = c(1, 1))

par(mfrow = c(3, 1))
hist(q3[1001:10000, 1], freq = F, main = "X1")
hist(q3[1001:10000, 2], freq = F, main = "X2")
hist(q3[1001:10000, 3], freq = F, main = "X3")
par(mfrow = c(1, 1))


# 4 ----

n = 100
mu_real = 5
tau_real = 1/4
y = rnorm(n, mu_real, sqrt(1/tau_real))

mu0 = 0
tau0 = 0.01
a = 0.01
b = 0.01

gibbs_normal = function(M, y, mu0, tau0, a, b) {
  n = length(y)
  ybar = mean(y)
  
  mu = NULL
  tau = NULL
  mu[1] = ybar
  tau[1] = 1 / var(y)
  
  for (t in 2:M) {
    # mu | tau, y
    prec = tau0 + n * tau[t-1]
    media = (tau0 * mu0 + n * tau[t-1] * ybar) / prec
    mu[t] = rnorm(1, media, sqrt(1 / prec))
    
    # tau | mu, y   (usa o mu que acabou de ser sorteado)
    tau[t] = rgamma(1, shape = a + n/2, rate = b + 0.5 * sum((y - mu[t])^2))
  }
  list(
    mu = mu,
    tau = tau
  )
}

q4 = gibbs_normal(10000, y, mu0, tau0, a, b)

mu_post = q4$mu[1001:10000]
tau_post = q4$tau[1001:10000]

par(mfrow = c(2, 1))
ts.plot(mu_post, ylab = expression(mu))
ts.plot(tau_post, ylab = expression(tau))
par(mfrow = c(1, 1))

par(mfrow = c(2, 1))
acf(mu_post)
acf(tau_post)
par(mfrow = c(1, 1))

par(mfrow = c(2, 1))
hist(mu_post, freq = F)
abline(v = mu_real, col = "red", lwd = 2)
hist(tau_post, freq = F)
abline(v = tau_real, col = "red", lwd = 2)
par(mfrow = c(1, 1))

plot(mu_post, tau_post, pch = 16, cex = 0.3, xlab = expression(mu), ylab = expression(tau))
