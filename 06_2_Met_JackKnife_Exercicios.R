# Aula 6.2 - Metodos de Jackknife - Exercicios

# 1 ----
x = c(2, 4, 6, 12)

n = length(x)

## 1.1 ----
x[-1]
x[-2]
x[-3]
x[-4]

## 1.2 ----
theta.i = NULL

for (i in 1:n) {
  
  theta.i[i] = mean(x[-i])
  
}

theta.i  # 7.333333 6.666667 6.000000 4.000000

## 1.3 ----
theta.bar = mean(theta.i)

theta.bar  # 6

## 1.4 ----

se.jack = sqrt(((n - 1)/n) * sum((theta.i - theta.bar)^2))

se.jack  # 2.160247

## 1.5 ----
se.trad = sd(x)/sqrt(n)

se.trad  # 2.160247

# 2 ----

x = c(1, 3, 5, 7, 9)

n = length(x)

## 2.1 ----

theta.hat = mean((x - mean(x))^2)

theta.hat  # 8

## 2.2 ----
theta.i = numeric(n)

for (i in 1:n) {
  
  x.sub = x[-i]
  
  theta.i[i] = mean((x.sub - mean(x.sub))^2)
  
}

theta.i  # 5.00 8.75 10.00 8.75 5.00

## 2.3 ----
theta.bar = mean(theta.i)

theta.bar  # 7.5

vies = (n - 1) * (theta.bar - theta.hat)

vies  # -2

## 2.4 ----
theta.jack = theta.hat - vies

theta.jack  # 10

var(x)  # 10

# 3 ----

## 1 ----

jackknife_cov = function(x, y){
  
  n = length(x)
  
  # Estimador da covariância
  estatistica = function(x, y){
    sum((x - mean(x)) * (y - mean(y)))/length(x)
  }
  
  theta_hat = estatistica(x, y)
  
  # Estimativas Jackknife (leave-one-out)
  theta_menosi = sapply(
    1:n,
    function(i){
      estatistica(x[-i], y[-i])
    }
  )
  
  theta_bar = mean(theta_menosi)
  
  # Estimador corrigido
  theta_jack = n*theta_hat - (n-1)*theta_bar
  
  # Vício estimado
  vies_jack = (n - 1) * (theta_bar - theta_hat)
  
  # Variância do estimador
  var_jack = ((n - 1) / n) * sum((theta_bar - theta_menosi)^2)
  
  # Erro padrão Jackknife
  se_jack = sqrt(var_jack)
  
  # Pseudovalores
  pseudovalores = n * theta_hat - (n - 1) * theta_menosi
  
  # Intervalo de confiança de 95%
  IC = c(
    theta_jack - qt(0.975, n - 1) * se_jack,
    theta_jack + qt(0.975, n - 1) * se_jack
  )
  
  list(
    estimativa_original = theta_hat,
    estimativa_jack = theta_jack,
    se_jack = se_jack,
    IC = IC
  )
}


## 2 ----

x = sort(sample(c(1:10), 10, replace = T))
y = sort(sample(c(1:10), 10, replace = T))

jackknife_cov(x, y)$se_jack

jackknife_cov(x, y)$IC

jackknife_cov(x, y)$estimativa_jack
cov(x, y)
