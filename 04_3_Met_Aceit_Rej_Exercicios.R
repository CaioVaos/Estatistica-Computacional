# Aula 4.3 - Aceitacao-Rejeicao - Exercicios

# 1 ----

## 1.1 ----

f <- function(x){
  4*x*I(x<=0.5)*I(x>0) + 4*(1-x)*I(x>0.5)*I(x<1)
}

curve(f)

k <- optimize(f, c(0,1), maximum = T)$objective

curve(dunif(x), add = T)
curve(k*dunif(x), add = T)

## 1.2 ----

rtriang.rej = function(n){
  
  x_gerado = NULL
  j = 0
  k = 2
  
  cont = 0
  while (j<n) {
    cont = cont +1
    
    y = runif(1)
    
    u = runif(1)
    
    px = f(y)
    qy = 1
    
    if(u<px/(k*qy)){
      j = j+1
      x_gerado[j] = y
    }
  }
  list(
    x_gerado = x_gerado,
    cont = cont
  )
}

set.seed(223054059)
x1 <- rtriang.rej(10000)
hist(x1$x_gerado, freq = F)
curve(f, add = T)

## 1.3 ----

# prop Teorica 
1/k

# prop Empirica
x1$cont #20186
10000/20186

## 1.4 ----

# prob empirica
mean(x1$x_gerado > 0.3 & x1$x_gerado < 0.7)

# prob teorica
integrate(f, 0.3, 0.7)$value

# 2 ----

## 2.1 ----

curve(dbeta(x,3,2))
curve(dnorm(x, 0.6, 0.35), add = T)

razao <- function(x){
  dbeta(x,3,2) / dnorm(x,0.6,0.35)
}

k <- optimize(razao, c(0,1),maximum = TRUE)$objective

curve(k*dnorm(x, 0.6, 0.35), add = T, col = 3)

## 2.2 ----

rbeta.rej = function(n){
  
  x_gerado = NULL
  y_todos  = NULL
  u_todos  = NULL
  aceito   = NULL
  j = 0
  
  cont = 0
  while (j<n) {
    cont = cont + 1
    
    y = rnorm(1, 0.6, 0.35)
    u = runif(1)
    
    px = dbeta(y, 3, 2)
    qy = dnorm(y, 0.6, 0.35)
    
    y_todos[cont] = y
    aceito[cont]  = (u < px/(k*qy))
    u_todos[cont] = u * k * qy
    
    if(aceito[cont]){
      j = j+1
      x_gerado[j] = y
    }
  }
  list(x_gerado = x_gerado, cont = cont,
       y_todos = y_todos, u_todos = u_todos, aceito = aceito)
}

set.seed(223054059)
x2 <- rbeta.rej(10000)
hist(x2$x_gerado, freq = F)
curve(dbeta(x, 3, 2), add = T, col = 2)

## 2.3 ----

plot(x2$y_todos, x2$u_todos,
     col = ifelse(x2$aceito, 3, 2),
     xlim = c(0,1))

mean(x2$x_gerado>0.8)

# 3 ----

## 3.1 ----

f <- function(x){
  ((2/pi)^(1/2))*exp(-((x^2)/2))*I(x>=0)
}

curve(f, xlim = c(0,4))
curve(dexp(x, 1), add = T, col = 2)

razao <- function(x){
  f(x)/dexp(x, 1)
}

k <- optimise(razao, c(0,1), maximum = T)$objective
curve(k*dexp(x, 1), add = T, col = 3)

## 3.2 ----

rhalfnormal.rej <- function(n){

  x_gerado = NULL
  j = 0
  cont = 0

  while(j<n){
    cont = cont + 1

    y = rexp(1, 1)

    px <- f(y)
    qy <- dexp(y)

    u = runif(1)

    if(u<(px/(k*qy))){
      j = j+1
      x_gerado[j] = y
    }
  }
  list(
    x_gerado = x_gerado,
    cont = cont
  )
}

set.seed(223054059)
x3 <- rhalfnormal.rej(10000)
hist(x3$x_gerado, probability = T)
curve(f, add = T, col = 2)

### 3.3 ----

# taxa teorica
1/k #0.7601672

# taxa empírica
length(x3$x_gerado)/x3$cont #0.7601735

### 3.4 ----

mean(x3$x_gerado)
var(x3$x_gerado)

