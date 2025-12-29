Model1_addcov<-function(n,p){
  rho = 0.2
  Sigma = matrix(rho,p-4,p-4)
  for(i in 1:(p-4)){
    Sigma[i,i] = 1
  }
  return(mvrnorm(n,rep(0,p-4),Sigma))
}

Model1_SRS<-function(n,alpha0,alpha1,betavec0,betavec1,pi,sigma0,sigma1,p){
  S0 = rep(0,n)
  S1 = rep(0,n)
  while(length(table(S0))!=4 || length(table(S1))!=4 || min(table(S0)) < 3 || min(table(S1)) < 3){
    X1 = rbeta(n,3,4)
    X2 = runif(n,-2,2)
    X3 = sample(c(1,-1),n,replace = TRUE, prob = c(0.5,0.5))
    X4 = sample(c(3,5),n,replace = TRUE, prob = c(0.6,0.4))
    S = sample(c(1,2,3,4),n,replace = TRUE, prob = c(0.2,0.3,0.3,0.2))
    A = SRS(n,pi)
    S0 = S[which(A == 0)]
    S1 = S[which(A == 1)]
  }
  Xbeta = cbind(X1,X2,X3,X4)
  Y0 = alpha0 + Xbeta%*%betavec0 + rnorm(n,sd = sigma0)
  Y1 = alpha1 + Xbeta%*%betavec1 + rnorm(n,sd = sigma1)
  X = cbind(X1,X2,X3,X4,Model1_addcov(n,p))
  Xrand = S
  return(list(A=A,S=S,X=data.frame(X),Xbeta = data.frame(Xbeta),Y=Y0*(1-A)+Y1*A))
}

Model1_WEI<-function(n,alpha0,alpha1,betavec0,betavec1,pi,sigma0,sigma1,p){
  S0 = rep(0,n)
  S1 = rep(0,n)
  while(length(table(S0))!=4 || length(table(S1))!=4 || min(table(S0)) < 3 || min(table(S1)) < 3){
    X1 = rbeta(n,3,4)
    X2 = runif(n,-2,2)
    X3 = sample(c(1,-1),n,replace = TRUE, prob = c(0.5,0.5))
    X4 = sample(c(3,5),n,replace = TRUE, prob = c(0.6,0.4))
    S = sample(c(1,2,3,4),n,replace = TRUE, prob = c(0.2,0.3,0.3,0.2))
    A = WEI(t(S),pi)
    S0 = S[which(A == 0)]
    S1 = S[which(A == 1)]
  }
  Xbeta = cbind(X1,X2,X3,X4)
  X = cbind(X1,X2,X3,X4,Model1_addcov(n,p))
  Y0 = alpha0 + Xbeta%*%betavec0 + rnorm(n,sd = sigma0)
  Y1 = alpha1 + Xbeta%*%betavec1 + rnorm(n,sd = sigma1)
  return(list(A=A,S=S,X=data.frame(X),Xbeta = data.frame(Xbeta),Y=Y0*(1-A)+Y1*A))
}

Model1_SBR<-function(n,alpha0,alpha1,betavec0,betavec1,pi,sigma0,sigma1,p){
  S0 = rep(0,n)
  S1 = rep(0,n)
  while(length(table(S0))!=4 || length(table(S1))!=4 || min(table(S0)) < 3 || min(table(S1)) < 3){
    X1 = rbeta(n,3,4)
    X2 = runif(n,-2,2)
    X3 = sample(c(1,-1),n,replace = TRUE, prob = c(0.5,0.5))
    X4 = sample(c(3,5),n,replace = TRUE, prob = c(0.6,0.4))
    S = sample(c(1,2,3,4),n,replace = TRUE, prob = c(0.2,0.3,0.3,0.2))
    A = SBR(t(S),pi)
    S0 = S[which(A == 0)]
    S1 = S[which(A == 1)]
  }
  Xbeta = cbind(X1,X2,X3,X4)
  X = cbind(X1,X2,X3,X4,Model1_addcov(n,p))
  Y0 = alpha0 + Xbeta%*%betavec0 + rnorm(n,sd = sigma0)
  Y1 = alpha1 + Xbeta%*%betavec1 + rnorm(n,sd = sigma1)
  return(list(A=A,S=S,X=data.frame(X),Xbeta = data.frame(Xbeta),Y=Y0*(1-A)+Y1*A))
}

Model1_BCD<-function(n,alpha0,alpha1,betavec0,betavec1,pi,sigma0,sigma1,p,lambda){
  S0 = rep(0,n)
  S1 = rep(0,n)
  while(length(table(S0))!=4 || length(table(S1))!=4 || min(table(S0)) < 3 || min(table(S1)) < 3){
    X1 = rbeta(n,3,4)
    X2 = runif(n,-2,2)
    X3 = sample(c(1,-1),n,replace = TRUE, prob = c(0.5,0.5))
    X4 = sample(c(3,5),n,replace = TRUE, prob = c(0.6,0.4))
    S = sample(c(1,2,3,4),n,replace = TRUE, prob = c(0.2,0.3,0.3,0.2))
    A = BCD(t(S),pi,lambda)
    S0 = S[which(A == 0)]
    S1 = S[which(A == 1)]
  }
  Xbeta = cbind(X1,X2,X3,X4)
  X = cbind(X1,X2,X3,X4,Model1_addcov(n,p))
  Y0 = alpha0 + Xbeta%*%betavec0 + rnorm(n,sd = sigma0)
  Y1 = alpha1 + Xbeta%*%betavec1 + rnorm(n,sd = sigma1)
  return(list(A=A,S=S,X=data.frame(X),Xbeta = data.frame(Xbeta),Y=Y0*(1-A)+Y1*A))
}

Model1_PS<-function(n,alpha0,alpha1,betavec0,betavec1,pi,sigma0,sigma1,p,weight,lambda){
  S0 = rep(0,n)
  S1 = rep(0,n)
  while(length(table(S0))!=4 || length(table(S1))!=4 || min(table(S0)) < 3 || min(table(S1)) < 3){
    X1 = rbeta(n,3,4)
    X2 = runif(n,-2,2)
    X3 = sample(c(1,-1),n,replace = TRUE, prob = c(0.5,0.5))
    X4 = sample(c(3,5),n,replace = TRUE, prob = c(0.6,0.4))
    S = sample(c(1,2,3,4),n,replace = TRUE, prob = c(0.2,0.3,0.3,0.2))
    A = PocSimue(t(S),pi,weight,lambda)
    S0 = S[which(A == 0)]
    S1 = S[which(A == 1)]
  }
  Xbeta = cbind(X1,X2,X3,X4)
  X = cbind(X1,X2,X3,X4,Model1_addcov(n,p))
  Y0 = alpha0 + Xbeta%*%betavec0 + rnorm(n,sd = sigma0)
  Y1 = alpha1 + Xbeta%*%betavec1 + rnorm(n,sd = sigma1)
  return(list(A=A,S=S,X=data.frame(X),Xbeta = data.frame(Xbeta),Y=Y0*(1-A)+Y1*A))
}

Model1_HH<-function(n,alpha0,alpha1,betavec0,betavec1,pi,sigma0,sigma1,p,omega,lambda){
  S0 = rep(0,n)
  S1 = rep(0,n)
  while(length(table(S0))!=4 || length(table(S1))!=4 || min(table(S0)) < 3 || min(table(S1)) < 3){
    X1 = rbeta(n,3,4)
    X2 = runif(n,-2,2)
    X3 = sample(c(1,-1),n,replace = TRUE, prob = c(0.5,0.5))
    X4 = sample(c(3,5),n,replace = TRUE, prob = c(0.6,0.4))
    S = sample(c(1,2,3,4),n,replace = TRUE, prob = c(0.2,0.3,0.3,0.2))
    A = HHue(t(S),pi,omega,lambda)
    S0 = S[which(A == 0)]
    S1 = S[which(A == 1)]
  }
  Xbeta = cbind(X1,X2,X3,X4)
  X = cbind(X1,X2,X3,X4,Model1_addcov(n,p))
  Y0 = alpha0 + Xbeta%*%betavec0 + rnorm(n,sd = sigma0)
  Y1 = alpha1 + Xbeta%*%betavec1 + rnorm(n,sd = sigma1)
  return(list(A=A,S=S,X=data.frame(X),Xbeta = data.frame(Xbeta),Y=Y0*(1-A)+Y1*A))
}
