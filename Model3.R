Model3_addcov<-function(n,p){
  rho = 0.5
  x = rep(1,p-2)
  for(i in 2:(p-2)){
    x[i] = x[i-1]*rho
  }
  Sigma = toeplitz(x)
  return(mvrnorm(n,rep(0,p-2),Sigma))
}

Model3_SRS<-function(n,alpha0,alpha1,betavec0,betavec1,pi,sigma0,sigma1,p){
  X1 = rbeta(n,3,4)
  X2 = runif(n,-2,2)
  X3 = rnorm(n,0,1)
  X4 = runif(n,0,2)
  S = sample(c(1,2),n,replace = TRUE, prob = c(0.4,0.6))
  A = SRS(n,pi)
  Xbeta = cbind(X1,X2)
  #Y0 = alpha0 + betavec0[1]*(X1*X2)/(X1+X2+2) + betavec0[2]*X1^2*(X2+X3) + rnorm(n,sd = sigma0)
  #Y1 = alpha1 + betavec1[1]*(X2+X4) + betavec1[2]*X2^2/exp(X1+2) + rnorm(n,sd = sigma1)
  Y0 = alpha0 + betavec0[1]*(X1*X2)/(X1+X2+2) + betavec0[2]*X1^2*(X2+X3) + sigma0 * rt(n,2)
  Y1 = alpha1 + betavec1[1]*(X2+X4) + betavec1[2]*X2^2/exp(X1+2) + sigma1 * rt(n,2)
  X = cbind(X1,X2,Model3_addcov(n,p))
  return(list(A=A,S=S,X=data.frame(X),Xbeta=data.frame(Xbeta),Y=Y0*(1-A)+Y1*A))
}

Model3_WEI<-function(n,alpha0,alpha1,betavec0,betavec1,pi,sigma0,sigma1,p){
  X1 = rbeta(n,3,4)
  X2 = runif(n,-2,2)
  X3 = rnorm(n,0,1)
  X4 = runif(n,0,2)
  S = sample(c(1,2),n,replace = TRUE, prob = c(0.4,0.6))
  A = WEI(t(S),pi)
  Xbeta = cbind(X1,X2)
  #Y0 = alpha0 + betavec0[1]*(X1*X2)/(X1+X2+2) + betavec0[2]*X1^2*(X2+X3) + rnorm(n,sd = sigma0)
  #Y1 = alpha1 + betavec1[1]*(X2+X4) + betavec1[2]*X2^2/exp(X1+2) + rnorm(n,sd = sigma1)
  Y0 = alpha0 + betavec0[1]*(X1*X2)/(X1+X2+2) + betavec0[2]*X1^2*(X2+X3) + sigma0 * rt(n,2)
  Y1 = alpha1 + betavec1[1]*(X2+X4) + betavec1[2]*X2^2/exp(X1+2) + sigma1 * rt(n,2)
  X = cbind(X1,X2,Model3_addcov(n,p))
  return(list(A=A,S=S,X=data.frame(X),Xbeta=data.frame(Xbeta),Y=Y0*(1-A)+Y1*A))
}

Model3_SBR<-function(n,alpha0,alpha1,betavec0,betavec1,pi,sigma0,sigma1,p){
  X1 = rbeta(n,3,4)
  X2 = runif(n,-2,2)
  X3 = rnorm(n,0,1)
  X4 = runif(n,0,2)
  S = sample(c(1,2),n,replace = TRUE, prob = c(0.4,0.6))
  A = SBR(t(S),pi)
  Xbeta = cbind(X1,X2)
  #Y0 = alpha0 + betavec0[1]*(X1*X2)/(X1+X2+2) + betavec0[2]*X1^2*(X2+X3) + rnorm(n,sd = sigma0)
  #Y1 = alpha1 + betavec1[1]*(X2+X4) + betavec1[2]*X2^2/exp(X1+2) + rnorm(n,sd = sigma1)
  Y0 = alpha0 + betavec0[1]*(X1*X2)/(X1+X2+2) + betavec0[2]*X1^2*(X2+X3) + sigma0 * rt(n,2)
  Y1 = alpha1 + betavec1[1]*(X2+X4) + betavec1[2]*X2^2/exp(X1+2) + sigma1 * rt(n,2)
  X = cbind(X1,X2,Model3_addcov(n,p))
  return(list(A=A,S=S,X=data.frame(X),Xbeta=data.frame(Xbeta),Y=Y0*(1-A)+Y1*A))
}

Model3_BCD<-function(n,alpha0,alpha1,betavec0,betavec1,pi,sigma0,sigma1,p,lambda){
  X1 = rbeta(n,3,4)
  X2 = runif(n,-2,2)
  X3 = rnorm(n,0,1)
  X4 = runif(n,0,2)
  S = sample(c(1,2),n,replace = TRUE, prob = c(0.4,0.6))
  A = BCD(t(S),pi,lambda)
  Xbeta = cbind(X1,X2)
  #Y0 = alpha0 + betavec0[1]*(X1*X2)/(X1+X2+2) + betavec0[2]*X1^2*(X2+X3) + rnorm(n,sd = sigma0)
  #Y1 = alpha1 + betavec1[1]*(X2+X4) + betavec1[2]*X2^2/exp(X1+2) + rnorm(n,sd = sigma1)
  Y0 = alpha0 + betavec0[1]*(X1*X2)/(X1+X2+2) + betavec0[2]*X1^2*(X2+X3) + sigma0 * rt(n,2)
  Y1 = alpha1 + betavec1[1]*(X2+X4) + betavec1[2]*X2^2/exp(X1+2) + sigma1 * rt(n,2)
  X = cbind(X1,X2,Model3_addcov(n,p))
  return(list(A=A,S=S,X=data.frame(X),Xbeta=data.frame(Xbeta),Y=Y0*(1-A)+Y1*A))
}

Model3_PS<-function(n,alpha0,alpha1,betavec0,betavec1,pi,sigma0,sigma1,p,weight,lambda){
  X1 = rbeta(n,3,4)
  X2 = runif(n,-2,2)
  X3 = rnorm(n,0,1)
  X4 = runif(n,0,2)
  S = sample(c(1,2),n,replace = TRUE, prob = c(0.4,0.6))
  A = PocSimue(t(S),pi,weight,lambda)
  Xbeta = cbind(X1,X2)
  #Y0 = alpha0 + betavec0[1]*(X1*X2)/(X1+X2+2) + betavec0[2]*X1^2*(X2+X3) + rnorm(n,sd = sigma0)
  #Y1 = alpha1 + betavec1[1]*(X2+X4) + betavec1[2]*X2^2/exp(X1+2) + rnorm(n,sd = sigma1)
  Y0 = alpha0 + betavec0[1]*(X1*X2)/(X1+X2+2) + betavec0[2]*X1^2*(X2+X3) + sigma0 * rt(n,2)
  Y1 = alpha1 + betavec1[1]*(X2+X4) + betavec1[2]*X2^2/exp(X1+2) + sigma1 * rt(n,2)
  X = cbind(X1,X2,Model3_addcov(n,p))
  return(list(A=A,S=S,X=data.frame(X),Xbeta=data.frame(Xbeta),Y=Y0*(1-A)+Y1*A))
}

Model3_HH<-function(n,alpha0,alpha1,betavec0,betavec1,pi,sigma0,sigma1,p,omega,lambda){
  X1 = rbeta(n,3,4)
  X2 = runif(n,-2,2)
  X3 = rnorm(n,0,1)
  X4 = runif(n,0,2)
  S = sample(c(1,2),n,replace = TRUE, prob = c(0.4,0.6))
  A = HHue(t(S),pi,omega,lambda)
  Xbeta = cbind(X1,X2)
  #Y0 = alpha0 + betavec0[1]*(X1*X2)/(X1+X2+2) + betavec0[2]*X1^2*(X2+X3) + rnorm(n,sd = sigma0)
  #Y1 = alpha1 + betavec1[1]*(X2+X4) + betavec1[2]*X2^2/exp(X1+2) + rnorm(n,sd = sigma1)
  Y0 = alpha0 + betavec0[1]*(X1*X2)/(X1+X2+2) + betavec0[2]*X1^2*(X2+X3) + sigma0 * rt(n,2)
  Y1 = alpha1 + betavec1[1]*(X2+X4) + betavec1[2]*X2^2/exp(X1+2) + sigma1 * rt(n,2)
  X = cbind(X1,X2,Model3_addcov(n,p))
  return(list(A=A,S=S,X=data.frame(X),Xbeta=data.frame(Xbeta),Y=Y0*(1-A)+Y1*A))
}