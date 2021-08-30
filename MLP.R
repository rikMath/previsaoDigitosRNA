rm(list=ls())

sech2<-function(u){
  return(((2/(exp(u) + exp(-u)))*(2/(exp(u) + exp(-u)))))
}

# Inicializa??o dos pesos
# Matriz Z nxp, neste cado n+1xp-> 3x2
Z<-matrix(runif(6)-0.5, ncol=2, nrow=3)

# Matriz W pxm, neste caso n+1xp-> 3x2
W<-matrix(runif(6)-0.5, ncol=2, nrow=3)

x<-matrix(c(0,0,0,1,1,0,1,1), ncol=2, byrow=T)
xatual<-matrix(nrow=3, ncol=1)
y<-matrix(c(-1,+1,+1,-1,+1,-1,-1,+1), ncol=2, byrow=T)

tol<-1.0
eta<-0.01
maxepocas<-2000
nepocas<-0
eepoca<-tol+1
N<-4
evec<-matrix(nrow=maxepocas, ncol=1)

while((nepocas < maxepocas) && (eepoca>tol)){
  ei2<-0
  # Sequ?ncia aleat?ria de treinamento
  xseq<-sample(N)
  for(i in 1:N){
    # Amostra dado da sequ?ncia aleat?ria
    irand<-xseq[i]
    xatual[(1:2), 1]<-x[irand,(1:2)]
    xatual[3,1]<-1

    yatual<-y[irand,(1:2)]

    U<-t(xatual)%*%Z # xatual ? 3x1 e Z ? 3x2
    H<-tanh(U)
    Haug<-cbind(H,1) # Haug ? 1x3

    O<-Haug%*%W
    yhat<-tanh(O)

    e<-yatual-yhat
    flinhaO<-sech2(O)
    dO<-e * flinhaO # Produto elemento a elemento

    Wminus<-W[-3,]  # Sa?da do vi?s n?o se propaga
    ehidden<-dO %*% t(Wminus) # dO 1x2 e W 3x2, ehidden ? 1x2
    flinhaU<-sech2(U)
    dU<-ehidden * flinhaU # Produto elemento a elemento

    W<-W+eta*(t(Haug) %*% dO)
    Z<-Z+eta*(xatual %*% dU)
    ei2<-ei2+(e %*% t(e))

  }
  # Incrementa n?mero de ?pocas
  nepocas<-nepocas+1
  evec[nepocas]<-ei2/N

  # Armazena erro por ?poca.
  eepoca<-evec[nepocas]
}
plot(evec[1:nepocas], type='l')
