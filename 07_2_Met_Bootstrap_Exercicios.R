# Aula 7.2 - Metodos de Bootstrap - Exercicios


# 1 ----

## 1 ----

# Prob de um elemento ser selecionado: 1/n
# Logo, de um elemento não ser selecionado: 1 - (1/n)
# em n amostras: [1 - (1/n)]^n

## 2 ----

# [1 - (1/n)]^n = [1 + (-1/n)]^n
# No limite:
# lim [1 + (-1/n)]^n = 1/e

## 3 ----

# Será o complementar da prob. anterior:
# 1 - (1/e)

# 2 ----

x_a <- c(12,15,18,22,25,28,33,40,45)

x_b <- c(8,10,11,12,14,15,17,18,20,22,25)


## 1 ----

# Estimativa pontual pelo metodo plug-in

theta_chapeu_plugin <- mean(x_a)/mean(x_b)

theta_chapeu_plugin


## 2 ----

bootstrap_razao = function(x_a, x_b, B){
  
  n_a = length(x_a)
  n_b = length(x_b)
  
  # estimativa pontual (baseada no estimador original)
  theta_hat = mean(x_a)/mean(x_b)
  
  # algoritmo Bootstrap
  replicas_theta_hat = numeric(B)
  
  for(b in 1:B) {
    
    reamostra_a = sample(x_a, size = n_a, replace = TRUE)
    
    reamostra_b = sample(x_b, size = n_b, replace = TRUE)
    
    replicas_theta_hat[b] =
      mean(reamostra_a)/mean(reamostra_b)
  }
  
  theta_bar = mean(replicas_theta_hat)
  
  # Vicio estimado
  vies_boot = theta_bar - theta_hat
  
  # Variancia do estimador
  var_boot = var(replicas_theta_hat)
  
  # Erro-padrao Bootstrap
  se_boot = sd(replicas_theta_hat)
  
  # Estimador corrigido
  theta_boot = theta_hat - vies_boot  
  
  list(
    estimativa_original = theta_hat,
    replicas_theta_hat = replicas_theta_hat,
    media_bootstrap = theta_bar,
    vies_boot = vies_boot,
    var_boot = var_boot,
    se_boot = se_boot,
    estimativa_boot = theta_boot
  )
}

bootstrap <- bootstrap_razao(x_a, x_b, 5000)

bootstrap$se_boot
bootstrap$vies_boot
bootstrap$estimativa_boot


## 3 ----

### Normal ----

bootstrap$estimativa_boot - qnorm(0.975) * bootstrap$se_boot
bootstrap$estimativa_boot + qnorm(0.975) * bootstrap$se_boot


### Empirico ----

quantile(bootstrap$replicas_theta_hat, 0.025)
quantile(bootstrap$replicas_theta_hat, 0.975)


# 3 ----

df <- mtcars

beta1 <- lm(mpg ~ hp, data = df)$coefficients[2]

beta1


bootstrap_ml = function(base, B){
  
  n = nrow(base)
  
  # estimativa pontual (baseada no estimador original)
  theta_hat =
    lm(mpg ~ hp, data = base)$coefficients[2]
  
  # algoritmo Bootstrap
  replicas_theta_hat = numeric(B)
  
  for(b in 1:B) {
    
    # Bootstrap de casos:
    # selecionamos n observacoes com reposicao
    
    indice = sample(1:n, size = n, replace = TRUE)
    
    reamostra = base[indice, ]
    
    # estimativa de beta1 na amostra Bootstrap
    
    replicas_theta_hat[b] =
      lm(mpg ~ hp, data = reamostra)$coefficients[2]
  }
  
  theta_bar = mean(replicas_theta_hat)
  
  # Vicio estimado
  vies_boot = theta_bar - theta_hat
  
  # Variancia do estimador
  var_boot = var(replicas_theta_hat)
  
  # Erro-padrao Bootstrap
  se_boot = sd(replicas_theta_hat)
  
  # Estimador corrigido
  theta_boot = theta_hat - vies_boot  
  
  # Intervalo de confianca Percentil
  IC_percentil =
    quantile(
      replicas_theta_hat,
      probs = c(0.025, 0.975)
    )
  
  list(
    estimativa_original = theta_hat,
    replicas_theta_hat = replicas_theta_hat,
    media_bootstrap = theta_bar,
    vies_boot = vies_boot,
    var_boot = var_boot,
    se_boot = se_boot,
    estimativa_boot = theta_boot,
    IC_percentil = IC_percentil
  )
}

bootstrap_ml_resultado <- bootstrap_ml(df, 2000)

bootstrap_ml_resultado$estimativa_original
bootstrap_ml_resultado$se_boot
bootstrap_ml_resultado$vies_boot
bootstrap_ml_resultado$estimativa_boot
bootstrap_ml_resultado$IC_percentil


# 4 ----

library(boot)

x <- c(12,15,18,20,25,30,32,38,42,45)
y <- c(100,115,110,130,140,150,168,175,190,200)

dados <- data.frame(
  x = x,
  y = y
)

cor_original <- cor(dados$x, dados$y)
cor_original

cor_bootstrap = function(data, indices){
  dados_boot = data[indices, ]
  cor(dados_boot$x, dados_boot$y)
}

bootstrap_cor = boot(
  data = dados,
  statistic = cor_bootstrap,
  R = 5000
)

bootstrap_cor$t0
bootstrap_cor$t
mean(bootstrap_cor$t)

vies_boot = mean(bootstrap_cor$t) - bootstrap_cor$t0
vies_boot

se_boot = sd(bootstrap_cor$t)
se_boot

IC_normal = boot.ci(
  bootstrap_cor,
  conf = 0.90,
  type = "norm"
)
IC_normal$normal

IC_basico = boot.ci(
  bootstrap_cor,
  conf = 0.90,
  type = "basic"
)
IC_basico$basic

IC_percentil = boot.ci(
  bootstrap_cor,
  conf = 0.90,
  type = "perc"
)
IC_percentil$percent
