# Aula 4.2 - Aceitacao-Rejeicao - Discretas

# 1 ----

f <- function(x){
  3*x^2
}

curve(f)
curve(dunif(x), col=2, add = T)

optimise(f, c(0,1), maximum = T) #2.999603
k <- 2.999603
curve(k*dunif(x), col=3, add = T)

r_quadr <- function(n){

  x_gerado = NULL
  j = 0
  k = 2.999603

  while(j<n){

    y = runif(1)

    u = runif(1)

    px = 3*y^2
    qy = 1

    if(u < px/(k*qy)){
      j = j+1
      x_gerado[j] = y
    }
  }
  x_gerado
}

x1 <- r_quadr(1000)
hist(x1, freq = F)
curve(f, add = T, col = 2)

# 2 ----

f <- function(x , a, b){
  dbeta(x , a, b)
}

a = 2
b = 4
curve(f(x, a, b))

optimize(f, c(0,1), a = a, b = b, maximum = T)$objective

curve(dunif(x))

r_beta <- function(n, a, b){
  x_gerado = NULL
  j = 0
  k = optimize(dbeta, c(0,1), shape1 = a, shape2 = b, maximum = T)$objective

  while(j<n){
    y = runif(1)
    u = runif(1)

    px = dbeta(y, a, b)
    qy = 1

    if(u < px/(k*qy)){
      j = j+1
      x_gerado[j] = y
    }
  }
  x_gerado
}

x2 <- r_beta(1000, 5, 20)
hist(x2,freq = F)
curve(dbeta(x, 5, 20), col = 2, add = T)

x2_R <- r_beta(1000, 5, 20)
hist(x2_R,freq = F, add = T, col= adjustcolor(4, alpha.f = 0.3))

summary(x2)
summary(x2_R)

boxplot(x2, x2_R)
