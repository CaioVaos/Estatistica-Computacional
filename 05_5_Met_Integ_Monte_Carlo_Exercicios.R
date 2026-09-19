# Aula 5.5 - Monte Carlo - Exercicios

# 1 ----

I.simples = function(n){
  x = runif(n)
  I = mean((1/(1+x^2))*exp(-x))
  return(I)
}

I.simples(10000)

simples <- NULL
for (i in 1:500) {
  simples[i] = I.simples(10000)
}

# 2 ----
curve(exp(-x)/(1+x^2), ylim = c(0,2))
curve(exp(x), add = T, col = 2)
curve(exp(-x), add = T, col = 3)

I.import = function(n){
  x = rexp(n)
  I = mean(I(x < 1)/(1 + x^2))
  return(I)
}

I.import(10000)

importancia <- NULL
for (i in 1:500) {
  importancia[i] = I.import(10000)
}

# 3 ----

I.ant = function(n){
  x = runif(n)
  y = (((1/(1+x^2))*exp(-x))+((1/(1+(1-x)^2))*exp(-(1-x))))/2
  I = mean(y)
  return(I)
}

I.ant(10000)

antitetica <- NULL
for (i in 1:500) {
  antitetica[i] = I.ant(10000)
}

# 4 ----

I.control = function(n){
  u = runif(n)
  x = (1/(1+u^2))*exp(-u)
  y = u
  
  x_cov = x - ((cov(x,y)/var(y)) * (y - 0.5))
  I = mean(x_cov)
  return(I)
}

I.control(10000)

controle <- NULL
for (i in 1:500) {
  controle[i] = I.control(10000)
}

# 5 ----
integrate(function(x) exp(-x)/(1+x^2), 0, 1)$value

# Comparação ----

# Comparação ----

var(simples)
var(importancia)
var(antitetica)
var(controle)

boxplot(list(simples, importancia, antitetica, controle))
