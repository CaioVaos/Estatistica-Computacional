# Aula 3.1 - Tranformação Inversa - Discreta

# Exercício 1 ----

runif.inv = function(n_esima, n = 1000){
  x = NULL

  for(i in 1:n){

    u = runif(1)
    prob_acumulada = 0

    for(j in 1:n_esima){

      prob_acumulada = prob_acumulada + 1/n_esima

      if(u <= prob_acumulada){
        x[i] = j
        break
      }
    }
  }
  x
}

x = runif.inv(10, n = 10000)
prop.table(table(x))

# Exercício 2 ----

rbin.inv <- function(n, p, N = 1000){
  
  x = c()

  for (i in 1:n) {
    u = runif(1)
    j = 0
    
    pj = (factorial(n)/(factorial(j)*factorial(n-j))) * p^j * (1-p)^(n-j)
    F = pj
    
    while(u >= F){
      j = j + 1
      pj = (factorial(n)/(factorial(j)*factorial(n-j))) * p^j * (1-p)^(n-j)
      F = F + pj
    }
    
    x[i] = j   
  }
  x
}

x = rbin.inv(100, 0.3)
prop.table(table(x))

mean(x)
100*0.3

# Exercício 3 ----
 
rgeom.inv = function(n, p){
  u = runif(n)
  x = floor(log(1-u)/log(1-p))
  x
}
 
rpois.inv = function(n, lambda){
  x = NULL
  for(i in 1:n){
    u = runif(1)
    j = 0
    p = exp(-lambda)
    F = p
    while(u >= F){
      j = j + 1
      p = (lambda/j) * p
      F = F + p
    }
    x[i] = j
  }
  x
}

# Uniforme
x1 = runif.inv(n = 1000, 6)
x1_R = runif(n = 1000, 0, 6)

par(mfrow=c(1,2))
hist(x1)
hist(x1_R)

# Geométrica
x2 = rgeom.inv(n = 1000, 0.3)
x2_R = rgeom(n = 1000, 0.3)

par(mfrow=c(1,2))
hist(x2)
hist(x2_R)

# Poisson
x3 = rpois.inv(n = 1000, 5)
x3_R = rpois(n = 1000, 5)

par(mfrow=c(1,2))
hist(x3)
hist(x3_R)

# Binomial
x4 = rbin.inv(N = 1000, 20, 0.3)
x4_R = rbinom(n = 1000, 20, 0.3)

par(mfrow=c(1,2))
hist(x4)
hist(x4_R)

