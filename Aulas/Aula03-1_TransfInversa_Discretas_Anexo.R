
#--------------------------------#
#---- Distribuição Geométrica ----
#--------------------------------#
rm(list=ls(all=TRUE))

n = 10000
p = 0.3

# Calculando a função de probabilidade teórica para os 10 primeiros valores da v.a..
# Lembrando que a parametrização do R utiliza o número de fracassos antes do 
# primeiro sucesso; enquanto utilizamos o número de realizações. 
# Deslocamos 1 em x para compatibilizar.

x = 0:9
f.p.teorica = dgeom(x,p)
f.p.teorica


#-----------------------------------------#
##---- Método da Transformação Inversa ----
#-----------------------------------------#

#-----------------------------#
###---- Forma mais simples ----
#-----------------------------#

rgeom.inv.simples = function(n,p){
  x = NULL  
  for (i in 1:n){
    u = runif(1)
    j = 1
    pj = p 
    F = pj
    while(u>=F){
      j = j+1
      pj = p*(1-p)^(j-1)
      F = F + pj
    }
    x[i] = j
  }
  x
}
  
x1 = rgeom.inv.simples(n,p)  
table(x1)[1:10]/n
f.p.teorica


#----------------------------------------------------------#
###---- Usando recursividade da função de probabilidade ----
#----------------------------------------------------------#

rgeom.inv.recurs = function(n,p){
  x = NULL  
  for (i in 1:n){
    u = runif(1)
    j = 1
    pj = p 
    F = pj
    while(u>=F){
      j = j+1
      pj = pj*(1-p)
      F = F + pj
    }
    x[i] = j
  }
  x
}

x2 = rgeom.inv.recurs(n,p)  
table(x2)[1:10]/n
f.p.teorica


#---------------------------------------#
###---- Usando transformação inversa ----
#---------------------------------------#

rgeom.inv.esperta = function(n,p){
  u = runif(n)
  x = floor(log(1-u)/log(1-p))
  x
}

x3 = rgeom.inv.esperta(n,p)
table(x3)[1:10]/n
f.p.teorica

#-----------------------------------------#
#---- Comparando as diferentes funções ----
#-----------------------------------------#

# Comparando as frequências dos valores gerados
round(cbind(f.p.teorica,table(x1)[1:10]/n,table(x2)[1:10]/n,table(x3)[1:10]/n),3)

# Comparando os tempos computacionais necessários para a geração

?system.time
?proc.time

n = 1000000
system.time(rgeom.inv.simples(n,p))
system.time(rgeom.inv.recurs(n,p))
system.time(rgeom.inv.esperta(n,p))
system.time(rgeom(n,p))


#--------------------------------#
#---- Distribuição Poisson ----
#--------------------------------#
rm(list=ls(all=TRUE))

n = 10000
lambda = 5

# Calculando a função de probabilidade teórica para os 10 primeiros valores da v.a..

x = 0:9
f.p.teorica = dpois(x,lambda)
f.p.teorica


#-----------------------------------------#
##---- Método da Transformação Inversa ----
#-----------------------------------------#

#-----------------------------#
###---- Forma mais simples ----
#-----------------------------#

rpois.inv.simples = function(n,lambda){
  x = NULL  
  for (i in 1:n){
    u = runif(1)
    j = 0
    pj = exp(-lambda) 
    F = pj
    while(u>=F){
      j = j+1
      pj = (lambda^j)*exp(-lambda)/factorial(j)
      F = F + pj
    }
    x[i] = j
  }
  x
}

x1 = rpois.inv.simples(n,lambda)  
table(x1)[1:10]/n
f.p.teorica


#----------------------------------------------------------#
###---- Usando recursividade da função de probabilidade ----
#----------------------------------------------------------#

rpois.inv.recurs = function(n,lambda){
  
  x = NULL  
  
  for (i in 1:n){
    u = runif(1)
    j = 0
    p = exp(-lambda)
    F = p
    while(u>=F){
      j = j+1
      p = (lambda/j)*p
      F = F + p
    }
    x[i] = j
  }
  x
}

x2 = rpois.inv.recurs(n,lambda)
table(x2)[1:10]/n
f.p.teorica

#-----------------------------------------#
#---- Comparando as diferentes funções ----
#-----------------------------------------#

# Comparando as frequências dos valores gerados
round(cbind(f.p.teorica,table(x1)[1:10]/n,table(x2)[1:10]/n),3)

# Comparando os tempos computacionais necessários para a geração

n = 1000000
system.time(rpois.inv.simples(n,lambda))
system.time(rpois.inv.recurs(n,lambda))
system.time(rpois(n,lambda))
