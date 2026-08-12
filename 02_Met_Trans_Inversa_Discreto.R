# Exercício 1 ----

rdunif.inv = function(n, nesima){
  x = NULL

  for(i in 1:n){

    u = runif(1)
    prob_acumulada = 0

    for(j in 1:nesima){

      prob_acumulada = prob_acumulada + 1/nesima

      if(u <= prob_acumulada){
        x[i] = j
        break
      }
    }
  }
  x
}

x = rdunif.inv(1000,5)
prop.table(table(x))

rdunif.inv = function(n,nesima){
  p = rep(1/nesima,nesima)
  F = cumsum(p)
  x = NULL
  for(i in 1:n){
    u = runif(1)
    x[i] = which(u<=F)[1]
  }
  x
}

x = rdunif.inv(1000,5)
prop.table(table(x))

# Exercício 2 ----

