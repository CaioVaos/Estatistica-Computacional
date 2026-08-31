# Aula 3.3 - Tranformação Inversa - Exercicios

# 1 ----

rexp.inv = function(n,lambda){
  u = runif(n)
  
  x = -log(1-u)/lambda
}

t1 = rexp.inv(10000,0.05)

hist(t1)

mean(t1)
mean(t1>30)
exp(-0.05*30)

# 2 ----

rpoli.inv = function(n){
  u = runif(n)
  
  x = (u)^(1/3)
}

poli = rpoli.inv(10000)

mean(poli)
var(poli)

# 3 ----
rpareto.inv = function(n,alpha,b){
  u = runif(n)
  
  x = b/(1-u)^(1/alpha)
}
 
sinistros = rpareto.inv(20000,2.5,1000)
 
hist(sinistros, breaks=100)
 
mean(sinistros)
quantile(sinistros,0.95)

# 4 ----
 
rnorm.boxmuller = function(n){
  u1 = runif(n)
  u2 = runif(n)
  
  z1 = sqrt(-2*log(u1))*cos(2*pi*u2)
  z2 = sqrt(-2*log(u1))*sin(2*pi*u2)
  
  z = c(z1,z2)
  z[1:n]
}
 
z = rnorm.boxmuller(10000)
 
hist(z, freq=FALSE)
curve(dnorm(x), add=TRUE, col="red")
 
mean(z>1.96)
1-pnorm(1.96)
 
system.time(rnorm(1e6))
system.time(rnorm.boxmuller(1e6))
 
rnorm.boxmuller2 = function(n,media=0,sd=1){
  media + sd*rnorm.boxmuller(n)
}
 
y = rnorm.boxmuller2(10000,5,2)
 
mean(y)
sd(y)