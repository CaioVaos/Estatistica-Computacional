# Aula 4.1 - Aceitacao-Rejeicao - Discretas

# 1 ----

rpois.rej = function(n,lambda, m = 1000){
  
  M = m
  k = (M+1)*max(dpois(0:M,lambda))
  
  x = rep(NA,n)
  
  cont = 0
  
  while(cont < n){
    
    # proposta uniforme
    y = sample(0:M,1,replace = TRUE)
    
    u = runif(1)
    px = dpois(y,lambda)
    gy = 1/(M+1)
    
    if(u < px/(k*gy)){
      cont = cont+1
      x[cont] = y
    }
    
  }
  x
}

x1 = rpois.rej(1000,7)

hist(x1, freq = FALSE)
points(0:max(x1), dpois(0:max(x1),7), col = "red", pch = 19)

mean(x1)
var(x1)


# 2 ----

rbinom.rej = function(n,size,prob){
  
  k = (size+1)*max(dbinom(0:size,size,prob))
  
  x = rep(NA,n)
  
  cont = 0
  
  while(cont < n){
    
    # proposta uniforme
    y = sample(0:size,1,replace = TRUE)
    
    u = runif(1)
    px = dbinom(y,size,prob)
    gy = 1/(size+1)
    
    if(u < px/(k*gy)){
      cont = cont+1
      x[cont] = y
    }
    
  }
  x
}

x2 = rbinom.rej(1000,20,0.8)

hist(x2, freq = FALSE)
points(0:max(x2), dbinom(0:max(x2),20, 0.8), col = "red", pch = 19)

mean(x2)
var(x2)
20*0.8
20*0.8*0.2
