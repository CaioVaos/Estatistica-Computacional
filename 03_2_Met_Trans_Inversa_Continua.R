# Aula 3.2 - Tranformação Inversa - Contínua

# Uniforme ----

runif.inv = function(n, a, b){
  u = runif(n)
  x = a + (b-a)*u
  x
}

x_runif = runif.inv(n, 3, 6)
x_runif_R = runif(n, min = 3, max = 6)

ks.test(x_runif, x_runif_R)

# Cauchy ----

# Cauchy ----

rcauchy.inv = function(n, mu, gamma){
  u = runif(n)
  x = mu + gamma*tan(pi*(u - 0.5))
  x
}

x_rcauchy = rcauchy.inv(n = 1000, mu = 2, gamma = 3)
x_rcauchy_R = rcauchy(n = 1000, location = 2, scale = 3)

ks.test(x_rcauchy, x_rcauchy_R)


# Pareto ----

rpareto.inv = function(n, a, b){
  u = runif(n)
  x = a*(1-u)^(-1/b)
  x
}

x_rpareto = rpareto.inv(n, a = 1, b = 2)
x_rpareto_R = VGAM::rpareto(n, scale = 1, shape = 2)

ks.test(x_rpareto, x_rpareto_R)

