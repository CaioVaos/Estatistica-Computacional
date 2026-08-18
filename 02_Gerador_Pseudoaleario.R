# Aula 2 - Geradores Congruenciais

# Setup ----

## Funções da aula

### x[n+1] = (a*x[n]) mod m
gcmulti = function(n, a, m, seed){
  x = numeric(n+1)
  x[1] = seed
  for(i in 2:(n+1))
    x[i] = (a*x[i-1]) %% m
  u = x/m
  u = u[-1]     # tira o x0, ficamos só com x1, x2, ..., xn
  u
}
 
# x[n+1] = (a*x[n]+c) mod m
gcmisto = function(n, a, c, m, seed){
  x = numeric(n+1)
  x[1] = seed
  for(i in 2:(n+1))
    x[i] = (a*x[i-1]+c) %% m
  u = x/m
  u = u[-1]
  u
}

# Exercício 1 ----

## 1) ----

u = gcmisto(40,7,5,32,1)
u
# > u
# [1] 0.37500 0.78125 0.62500 0.53125 0.87500 0.28125 0.12500 0.03125 0.37500 0.78125
# [11] 0.62500 0.53125 0.87500 0.28125 0.12500 0.03125 0.37500 0.78125 0.62500 0.53125
# O pérido é 16: a cada 16 números a seguência torna a se repitir

## 2) ----

u1 = gcmisto(1000,2,3,15,123)
u2 = gcmisto(5000,1103515245,12345,2^31,123)

hist(u1)
hist(u2)

# o u1 fica com barras esparssas 
# enquato que u2 fica mais proximo a uma U(0,1)

## 3) ----

u = gcmisto(5000,1103515245,12345,2^31,123)
mean(u)

# > mean(u)
# [1] 0.4951754
# a media da pratica fica bem proxima da media da teorica,
# Sinal de um pro gerador

## 4) ----

u = gcmisto(5000,1103515245,12345,2^31,123)
plot(u[-5000],u[-1],
     pch=19,
     cex=.2)
cor(u[-1],u[-length(u)])

# [1] -0.01912896
# Os pontos preenchem bem o quadrante de fomra aleatoria, sinal de independência
# Correlação proxíma do 0, sinal de independência

## 5) ----

?RNGkind
RNGkind()

# O padrão do R desde a versão 1.7.0 é o gerador "Mersenne-Twister",
# criado por Matsumoto e Nishimura (1998). Ele tem um período
# astronomicamente maior (2^19937 - 1) e propriedades estatísticas
# muito melhores do que um gerador congruencial simples

# Exercício 2 ----

## 1) ----

# O gerador de Park-Miller é um Gerador Congruencial MULTIPLICATIVO
# com a = 7^5 = 16807 e m = 2^31 - 1

parkmiller = function(n,seed){
  gcmulti(n,7^5,2^31-1,seed)
}

u = parkmiller(10000,12345)
hist(u)

## 2) ----

u_gcl = gcmulti(10000,16807,2^31-1,12345)
u_gcm = gcmisto(10000,16807,2,2^31-1,12345)

par(mfrow=c(1,2))
hist(u_gcl, main="GCL (multiplicativo)")
hist(u_gcm, main="GCM (misto)")

# Com bons parâmetros, os dois tipos de gerador produzem
# histogramas parecidos com o de uma distribuição uniforme

## 3) ----

# Autocorrelação de ordem k mede a correlação entre u[i] e u[i+k].
# Se os números forem aproximadamente independentes, essas
# correlações devem ficar próximas de zero.

u = gcmisto(10000,16807,2,2^31-1,12345)
cor(u[-1],u[-length(u)])
cor(u[-2],u[-length(u)])
cor(u[-5],u[-length(u)])
cor(u[-10],u[-length(u)])

u = gcmisto(40,7,5,32,1)

cor(u[-1],u[-length(u)])
cor(u[-2],u[-length(u)])
cor(u[-5],u[-length(u)])
cor(u[-10],u[-length(u)])
for (i in 10:30) {
  print(paste(i, cor(u[-i],u[-length(u)])))
}


## 4) ----

breaks = seq(0,1,length=11)

for(n in c(100,1000,10000,100000)){
  u = gcmisto(n,16807,2,2^31-1,12345)
  obs = table(cut(u,breaks=breaks))
  print(chisq.test(obs))
}

## 5) -----
# Além do teste qui-quadrado de aderência, existem outros testes
# clássicos usados para avaliar geradores de números
# pseudoaleatórios, entre eles:
#
# - Teste de Kolmogorov-Smirnov (compara a f.d.a. empírica com a
#   f.d.a. teórica da Uniforme(0,1));
# - Teste dos Runs / Teste de Sequências (verifica padrões de
#   subida e descida na sequência gerada);
# - Teste Espectral (avalia a estrutura de rede formada por pares
#   ou triplas sucessivas de números gerados);
# - Teste do Poker (agrupa os dígitos dos números gerados e
#   verifica se certas combinações aparecem com a frequência
#   esperada);
# - Bateria de testes "Diehard" e "TestU01" (conjuntos completos
#   de testes estatísticos usados para avaliar geradores
#   modernos, incluindo o Mersenne-Twister).
