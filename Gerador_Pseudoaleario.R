# Exercício 1 ----

## 1) ----

gcmisto = function(n,a,c,m,seed){

  x = numeric(n+1)
  x[1] = seed

  for(i in 2:(n+1))
    x[i] = (a*x[i-1]+c) %% m

  u = x/m
  u = u[-1]
}

u = gcmisto(40,7,5,32,1)
u

## 2) ----

u1 = gcmisto(1000,2,3,15,123)
u2 = gcmisto(5000,1103515245,12345,2^31,123)

hist(u1)
hist(u2)

## 3) ----

u = gcmisto(5000,1103515245,12345,2^31,123)
mean(u)

## 4) ----

plot(u[-length(u)],u[-1],
     pch=19,
     cex=.2)

cor(u[-1],u[-length(u)])



## 5) ----

?RNGkind

# Exercício 2 ----

gcmulti = function(n,a,m,seed){

  x = numeric(n+1)
  x[1] = seed

  for(i in 2:(n+1))
    x[i] = (a*x[i-1]) %% m

  u = x/m
  u = u[-1]
}

gcmisto = function(n,a,c,m,seed){

  x = numeric(n+1)
  x[1] = seed

  for(i in 2:(n+1))
    x[i] = (a*x[i-1]+c) %% m

  u = x/m
  u = u[-1]
}

## 1) ----

parkmiller = function(n,seed){
  gcmulti(n,16807,2^31-1,seed)
}

u = parkmiller(10000,12345)
hist(u)

## 2) ----

u_gcl = gcmulti(10000,16807,2^31-1,12345)
u_gcm = gcmisto(10000,16807,2,2^31-1,12345)

par(mfrow=c(1,2))
hist(u_gcl, main="GCL (multiplicativo)")
hist(u_gcm, main="GCM (misto)")

## 3) ----

u = gcmisto(10000,16807,2,2^31-1,12345)
u
cor(u[-1],u[-length(u)])
cor(u[-2],u[-length(u)])
cor(u[-5],u[-length(u)])
cor(u[-10],u[-length(u)])

u = gcmisto(40,7,5,32,1)
u
cor(u[-1],u[-length(u)])
cor(u[-2],u[-length(u)])
cor(u[-5],u[-length(u)])
cor(u[-10],u[-length(u)])
cor(u[-17],u[-length(u)])
cor(u[-18],u[-length(u)])
cor(u[-19],u[-length(u)])

## 4) ----

breaks = seq(0,1,length=11)

for(n in c(100,1000,10000,100000)){
  u = gcmisto(n,16807,2,2^31-1,12345)
  obs = table(cut(u,breaks=breaks))
  print(chisq.test(obs))
}

## 5) -----
