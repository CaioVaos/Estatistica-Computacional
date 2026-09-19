# Aula 5.3 - Monte Carlo - Exercicios

# 1 ----

circulo.MC <- function(n){
  x = runif(n, -1, 1)
  
  g = (1-x^2)^(1/2)
  
  2*mean(g)
}

circulo.MC(1000)
pi/2

Ns = c(100, 1000, 10000, 100000)
estimativas = sapply(Ns, circulo.MC)
erros = abs(estimativas - pi/2)

data.frame(N = Ns, estimativa = estimativas, erro = erros)

# 2 ----

exp.MC <- function(n){
   x = runif(n, 0, 1)
   
   g = exp(-x^2)/1
   
   mean(g)
}

exp.MC(10000)

integrate(function(x) exp(-x^2), 0, 1)

x <- seq(0, 2, by = 0.01)
f <- exp(-x^2)
plot(x, f, type = "l", xlim = c(0,2), ylim = c(0,1))
x2 <- seq(0, 1, by = 0.01)
f2 <- exp(-x2^2)
y2 <- rep(0, length(x2))
polygon(
  c(x2, rev(x2)), 
  c(f2, rev(y2)),
  col = "lightblue",
  border = NA
)

# 3 ----
exp.bi.MC <- function(n) {
  x <- runif(n, 0, 1)
  y <- runif(n, 0, 1)
  
  g <- exp(-(x^2 + y^2))
  
  mean(g)
}

exp.bi.MC(10000)

a <- integrate(function(x) exp(-x^2), 0, 1)$value
b <- integrate(function(x) exp(-x^2), 0, 1)$value
a*b

# 4 ----
esfera.MC <- function(n) {
  x = runif(n, -1, 1)
  y = runif(n, -1, 1)
  z = runif(n, -1, 1)
  
  indicadora = (x^2 + y^2 + z^2 <= 1)
  
  list(
    integral = 8*mean(indicadora),
    pontos = indicadora
  )
}

esfera.MC(10000)$integral

(4*pi)/3

# 5 ----
gama.MC <-  function(n,lambda){
  x = rexp(n, lambda)
  
  g = ((x^2)*(exp(-x)))/dexp(x, lambda)
  
  mean(g)
}

integrate(function(x) (x^2)*exp(-x), 0, Inf)$value

gama.MC(10000, 0.5)
gama.MC(10000, 1)
gama.MC(10000, 2)

gama.MC_0.5 <- NULL
gama.MC_1 <- NULL
gama.MC_2 <- NULL
for (i in 1:500){
  gama.MC_0.5[i] = gama.MC(10000, 0.5)
  gama.MC_1[i] = gama.MC(10000, 0.5)
  gama.MC_2[i] = gama.MC(10000, 0.5)
}
data.frame(
  lambda     = c(0.5, 1, 2),
  estimativa = c(
    mean(gama.MC_0.5),
    mean(gama.MC_1),
    mean(gama.MC_2)
  ),
  variancia  = c(
    var(gama.MC_0.5),
    var(gama.MC_1),
    var(gama.MC_2)
  )
)
boxplot(list(gama.MC_0.5, gama.MC_1, gama.MC_2))
