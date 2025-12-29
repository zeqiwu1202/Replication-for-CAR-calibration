library(mda)
library(splines)
library(randomForest)
library(caret)
library(MASS)
library(glmnet)
library(kernlab)
library(nnet)
library(gbm)
library(rpart)
library(caretEnsemble)
library(elasticnet)
library(np)
library(earth)
library(neuralnet)

#Stratification
Strat<-function(X3,X4){
  S = rep(0,length(X3))
  S[which(X3 == 1 & X4 == 3)] = 1
  S[which(X3 == 1 & X4 == 5)] = 2
  S[which(X3 == -1 & X4 == 3)] = 3
  S[which(X3 == -1 & X4 == 5)] = 4
  return(S)
}

stratum_cent<-function(A,S,X,Y,strt,strt_num){
  for(i in 1:strt_num){
    indk = which(S == strt[i])
    X[indk,] = X[indk,] - matrix(rep(apply(X[indk,],1,mean),n),ncol = n)
    Y[indk] = Y[indk] - mean(Y[indk])
  }
  return(data.frame(data.Y = Y,X))
}

samplevar<-function(A,S,R,strt){
  strt_num = length(strt)
  pait = length(which(A == 1))/length(A)
  r0 = rep(0,strt_num)
  r1 = rep(0,strt_num)
  nk0 = rep(0,strt_num)
  nk1 = rep(0,strt_num)
  sr0 = rep(0,strt_num)
  sr1 = rep(0,strt_num)
  for(i in 1:strt_num){
    inds = which(S == strt[i])
    nk1[i] = sum(A[inds])
    nk0[i] = sum(1 - A[inds])
    r0[i] = sum(R[inds]*(1-A[inds]))/nk0[i]
    r1[i] = sum(R[inds]*A[inds])/nk1[i]
    sr0[i] = sum((R[inds]-r0[i])^2*(1-A[inds]))/nk0[i]
    sr1[i] = sum((R[inds]-r1[i])^2*A[inds])/nk1[i]
  }
  pn = (nk0+nk1)/sum(nk0+nk1)
  r1 = r1 - sum(A*R)/sum(nk1)
  r0 = r0 - sum((1-A)*R)/sum(nk0)
  return(sum(pn*sr1)/pait+sum(pn*sr0)/(1-pait)+sum(pn*(r1-r0)^2))
}

est<-function(A,S,R,R_strat,strt,strt_num,n){
  pn = numeric(strt_num)
  Yk = numeric(strt_num)
  Yk_strat = numeric(strt_num)
  for(i in 1:strt_num){
    indk1 = which(A == 1 & S == strt[i])
    indk0 = which(A == 0 & S == strt[i])
    Yk_strat[i] = mean(R_strat[indk1]) - mean(R_strat[indk0])
    Yk[i] = mean(R[indk1]) - mean(R[indk0])
    pn[i] = (length(indk1)+length(indk0))/n
  }
  tau_est = sum(pn*Yk)
  var_est = samplevar(A,S,R,strt)
  tau_strat = sum(pn*Yk_strat)
  var_strat = samplevar(A,S,R_strat,strt)
  res = matrix(c(tau_est,sqrt(var_est/n),tau_strat,sqrt(var_strat/n)),nrow = 2,ncol = 2,byrow = TRUE)
  return(res)
}

est_lin<-function(data){
  n = length(data$A)
  strt = unique(data$S)
  strt_num = length(strt)
  pai = numeric(strt_num)
  ind0 = which(data$A == 0)
  ind1 = which(data$A == 1)
  datafit = data.frame(data$Y,data$Xbeta)
  fit1 = lm(data.Y~.,data = datafit[ind1,])
  fit0 = lm(data.Y~.,data = datafit[ind0,])
  h1 = predict(fit1,data$Xbeta)
  h0 = predict(fit0,data$Xbeta)
  pait = length(which(data$A == 1))/n
  R = data$Y - ((1-pait)*h1+pait*h0)
  R_strat = data$Y
  for(i in 1:strt_num){
    indk1 = which(data$A == 1 & data$S == strt[i])
    indk0 = which(data$A == 0 & data$S == strt[i])
    fitk1 = lm(data.Y~.,data = datafit[indk1,])
    fitk0 = lm(data.Y~.,data = datafit[indk0,])
    indk = c(indk1,indk0)
    hk1 = predict(fitk1,datafit[indk,])
    hk0 = predict(fitk0,datafit[indk,])
    pai[i] = length(indk1)/length(indk)
    R_strat[indk] = R_strat[indk] - ((1-pai[i])*hk1+pai[i]*hk0)
  }
  res = est(data$A,data$S,R,R_strat,strt,strt_num,n)
  return(res)
}

est_lin_cent<-function(data){
  n = length(data$A)
  strt = unique(data$S)
  strt_num = length(strt)
  pai = numeric(strt_num)
  ind0 = which(data$A == 0)
  ind1 = which(data$A == 1)
  datafit = stratum_cent(data$A,data$S,data$Xbeta,data$Y,strt,strt_num)
  fit1 = lm(data.Y~.,data = datafit[ind1,])
  fit0 = lm(data.Y~.,data = datafit[ind0,])
  h1 = predict(fit1,datafit[,-1])
  h0 = predict(fit0,datafit[,-1])
  pait = length(which(data$A == 1))/n
  R = data$Y - ((1-pait)*h1+pait*h0)
  R_strat = data$Y
  for(i in 1:strt_num){
    indk1 = which(data$A == 1 & data$S == strt[i])
    indk0 = which(data$A == 0 & data$S == strt[i])
    fitk1 = lm(data.Y~.,data = datafit[indk1,])
    fitk0 = lm(data.Y~.,data = datafit[indk0,])
    indk = c(indk1,indk0)
    hk1 = predict(fitk1,datafit[indk,-1])
    hk0 = predict(fitk0,datafit[indk,-1])
    pai[i] = length(indk1)/length(indk)
    R_strat[indk] = R_strat[indk] - ((1-pai[i])*hk1+pai[i]*hk0)
  }
  res = est(data$A,data$S,R,R_strat,strt,strt_num,n)
  return(res)
}

kernel_process<-function(Xbeta){
  for(i in 1:ncol(Xbeta)){
    if(length(table(Xbeta[,i]))<=10){
      Xbeta[,i] <- factor(Xbeta[,i])
    }
  }
  return(Xbeta)
}

est_kernel<-function(data){
  n = length(data$A)
  strt = unique(data$S)
  strt_num = length(strt)
  pai = numeric(strt_num)
  ind0 = which(data$A == 0)
  ind1 = which(data$A == 1)
  Xbeta = kernel_process(data$Xbeta)
  Xp = ncol(Xbeta)
  bw1 <- npregbw(xdat = Xbeta[ind1,],ydat = data$Y[ind1],regtype = "ll",bwmethod="cv.aic")
  bw0 <- npregbw(xdat = Xbeta[ind0,],ydat = data$Y[ind0],regtype = "ll",bwmethod="cv.aic")
  fit1 = npreg(txdat = Xbeta[ind1,],tydat = data$Y[ind1],exdat = Xbeta,regtype = "ll",bws = bw1)
  fit0 = npreg(txdat = Xbeta[ind0,],tydat = data$Y[ind0],exdat = Xbeta,regtype = "ll",bws = bw0)
  h1 = fitted(fit1)
  h0 = fitted(fit0)
  pait = length(which(data$A == 1))/n
  R = data$Y - ((1-pait)*h1+pait*h0)
  R_strat = data$Y
  for(i in 1:strt_num){
    indk1 = which(data$A == 1 & data$S == strt[i])
    indk0 = which(data$A == 0 & data$S == strt[i])
    indk = c(indk1,indk0)
    bwk1 <- npregbw(xdat = Xbeta[indk1,],ydat = data$Y[indk1],regtype = "ll",bwmethod="cv.aic")
    bwk0 <- npregbw(xdat = Xbeta[indk0,],ydat = data$Y[indk0],regtype = "ll",bwmethod="cv.aic")
    fitk1 = npreg(txdat = Xbeta[indk1,],tydat = data$Y[indk1],exdat = Xbeta[indk,],regtype = "ll",bws = bwk1)
    fitk0 = npreg(txdat = Xbeta[indk0,],tydat = data$Y[indk0],exdat = Xbeta[indk,],regtype = "ll",bws = bwk0)
    hk1 = fitted(fitk1)
    hk0 = fitted(fitk0)
    pai[i] = length(indk1)/length(indk)
    R_strat[indk] = R_strat[indk] - ((1-pai[i])*hk1+pai[i]*hk0)
  }
  res = est(data$A,data$S,R,R_strat,strt,strt_num,n)
  return(res)
}

est_ns<-function(data,Model){
  n = length(data$A)
  strt = unique(data$S)
  strt_num = length(strt)
  pai = numeric(strt_num)
  ind0 = which(data$A == 0)
  ind1 = which(data$A == 1)
  datafit = data.frame(data$Y,data$Xbeta)
  if(Model == "Model1"||Model == "Model5"){
    fit1 = lm(data.Y~ns(X1,10)+ns(X2,10)+X3+X4,
                    data = datafit[ind1,])
    fit0 = lm(data.Y~ns(X1,10)+ns(X2,10)+X3+X4,
                    data = datafit[ind0,])
    h0 = predict(fit0,data$Xbeta)
    h1 = predict(fit1,data$Xbeta)
    pait = length(which(data$A == 1))/n
    R = data$Y - ((1-pait)*h1+pait*h0)
    R_strat = data$Y
    for(i in 1:strt_num){
      indk1 = which(data$A == 1 & data$S == strt[i])
      indk0 = which(data$A == 0 & data$S == strt[i])
      fitk1 = lm(data.Y~ns(X1,10)+ns(X2,10)+X3+X4,
                data = datafit[indk1,])
      fitk0 = lm(data.Y~ns(X1,10)+ns(X2,10)+X3+X4,
                data = datafit[indk0,])
      indk = c(indk1,indk0)
      hk1 = predict(fitk1,data$Xbeta[indk,])
      hk0 = predict(fitk0,data$Xbeta[indk,])
      pai[i] = length(indk1)/length(indk)
      R_strat[indk] = R_strat[indk] - ((1-pai[i])*hk1+pai[i]*hk0)
    }
    res = est(data$A,data$S,R,R_strat,strt,strt_num,n)
    return(res)
  }
  else if(Model == "Model3"||Model == "Model7"){
    fit1 = lm(data.Y~ns(X1,10)+ns(X2,10)+ns(X3,10)+ns(X4,10),
              data = datafit[ind1,])
    fit0 = lm(data.Y~ns(X1,10)+ns(X2,10)+ns(X3,10)+ns(X4,10),
              data = datafit[ind0,])
    h0 = predict(fit0,data$Xbeta)
    h1 = predict(fit1,data$Xbeta)
    pait = length(which(data$A == 1))/n
    R = data$Y - ((1-pait)*h1+pait*h0)
    R_strat = data$Y
    for(i in 1:strt_num){
      indk1 = which(data$A == 1 & data$S == strt[i])
      indk0 = which(data$A == 0 & data$S == strt[i])
      fitk1 = lm(data.Y~ns(X1,10)+ns(X2,10)+X3+X4,
                 data = datafit[indk1,])
      fitk0 = lm(data.Y~ns(X1,10)+ns(X2,10)+X3+X4,
                 data = datafit[indk0,])
      indk = c(indk1,indk0)
      hk1 = predict(fitk1,data$Xbeta[indk,])
      hk0 = predict(fitk0,data$Xbeta[indk,])
      pai[i] = length(indk1)/length(indk)
      R_strat[indk] = R_strat[indk] - ((1-pai[i])*hk1+pai[i]*hk0)
    }
    res = est(data$A,data$S,R,R_strat,strt,strt_num,n)
    return(res)
  }
  else if(Model == "Model2"||Model == "Model4"||Model == "Model6"||Model == "Model8"){
    fit1 = lm(data.Y~ns(X1,10)+ns(X2,10),
                    data = datafit[ind1,])
    fit0 = lm(data.Y~ns(X1,10)+ns(X2,10),
                    data = datafit[ind0,])
    h0 = predict(fit0,data$Xbeta)
    h1 = predict(fit1,data$Xbeta)
    pait = length(which(data$A == 1))/n
    R = data$Y - ((1-pait)*h1+pait*h0)
    R_strat = data$Y
    for(i in 1:strt_num){
      indk1 = which(data$A == 1 & data$S == strt[i])
      indk0 = which(data$A == 0 & data$S == strt[i])
      fitk1 = lm(data.Y~ns(X1,10)+ns(X2,10),
                data = datafit[indk1,])
      fitk0 = lm(data.Y~ns(X1,10)+ns(X2,10),
                data = datafit[indk0,])
      indk = c(indk1,indk0)
      hk1 = predict(fitk1,data$Xbeta[indk,])
      hk0 = predict(fitk0,data$Xbeta[indk,])
      pai[i] = length(indk1)/length(indk)
      R_strat[indk] = R_strat[indk] - ((1-pai[i])*hk1+pai[i]*hk0)
    }
    res = est(data$A,data$S,R,R_strat,strt,strt_num,n)
    return(res)
  }
}

est_MARS<-function(data){
  n = length(data$A)
  strt = unique(data$S)
  strt_num = length(strt)
  pai = numeric(strt_num)
  ind0 = which(data$A == 0)
  ind1 = which(data$A == 1)
  datafit = data.frame(data$Y,data$Xbeta)
  fit1 = earth(data.Y~.,data = datafit[ind1,])
  fit0 = earth(data.Y~.,data = datafit[ind0,])
  h1 = predict(fit1,data$Xbeta)
  h0 = predict(fit0,data$Xbeta)
  pait = length(ind1)/n
  R = data$Y - ((1-pait)*h1+pait*h0)
  R_strat = data$Y
  for(i in 1:strt_num){
    indk1 = which(data$A == 1 & data$S == strt[i])
    indk0 = which(data$A == 0 & data$S == strt[i])
    fitk1 = earth(data.Y~.,data = datafit[indk1,])
    fitk0 = earth(data.Y~.,data = datafit[indk0,])
    indk = c(indk1,indk0)
    hk1 = predict(fitk1,datafit[indk,])
    hk0 = predict(fitk0,datafit[indk,])
    pai[i] = length(indk1)/length(indk)
    R_strat[indk] = R_strat[indk] - ((1-pai[i])*hk1+pai[i]*hk0)
  }
  res = est(data$A,data$S,R,R_strat,strt,strt_num,n)
  return(res)
}

#high-dimensional methods
est_lasso<-function(data){
  n = length(data$A)
  strt = unique(data$S)
  strt_num = length(strt)
  pai = numeric(strt_num)
  ind0 = which(data$A == 0)
  ind1 = which(data$A == 1)
  data$X = as.matrix(data$X)
  fit0<-cv.glmnet(data$X[ind0,],data$Y[ind0],family = "gaussian",standardize = FALSE,nfolds = 5)
  fit1<-cv.glmnet(data$X[ind1,],data$Y[ind1],family = "gaussian",standardize = FALSE,nfolds = 5)
  beta0<-coef(fit0,s = fit0$lambda.1se)[-1]
  beta1<-coef(fit1,s = fit1$lambda.1se)[-1]
  pait = length(which(data$A == 1))/n
  R = data$Y - data$X%*%((1-pait)*beta1+pait*beta0)
  R_strat = data$Y
  for(i in 1:strt_num){
    indk1 = which(data$A == 1 & data$S == strt[i])
    indk0 = which(data$A == 0 & data$S == strt[i])
    fitk0<-cv.glmnet(data$X[indk0,],data$Y[indk0],family = "gaussian",standardize = FALSE,nfolds = 5)
    fitk1<-cv.glmnet(data$X[indk1,],data$Y[indk1],family = "gaussian",standardize = FALSE,nfolds = 5)
    betak0<-coef(fitk0,s = fitk0$lambda.1se)[-1]
    betak1<-coef(fitk1,s = fitk1$lambda.1se)[-1]
    indk = c(indk1,indk0)
    pai[i] = length(indk1)/length(indk)
    R_strat[indk] = R_strat[indk] - data$X[indk,]%*%((1-pai[i])*betak1+pai[i]*betak0)
  }
  res = est(data$A,data$S,R,R_strat,strt,strt_num,n)
  return(res)
}

est_nnet<-function(data){
  n = length(data$A)
  strt = unique(data$S)
  strt_num = length(strt)
  pai = numeric(strt_num)
  ind0 = which(data$A == 0)
  ind1 = which(data$A == 1)
  datafit = data.frame(data$Y,data$X)
  Y.range = max(data$Y) - min(data$Y)
  Y.min = min(data$Y)
  datafit$data.Y <- (datafit$data.Y - Y.min)/Y.range
  fit1 = nnet(data.Y~.,data = datafit[ind1,],size = 2,linout = TRUE,rang = 0.1,decay = 0.02, maxit = 1000,maxNwt = 1500)
  fit0 = nnet(data.Y~.,data = datafit[ind0,],size = 2,linout = TRUE,rang = 0.1,decay = 0.02, maxit = 1000,maxNwt = 1500)
  h1 = predict(fit1,data$X)*Y.range + Y.min
  h0 = predict(fit0,data$X)*Y.range + Y.min
  pait = length(which(data$A == 1))/n
  R = data$Y - ((1-pait)*h1+pait*h0)
  R_strat = data$Y
  for(i in 1:strt_num){
    indk1 = which(data$A == 1 & data$S == strt[i])
    indk0 = which(data$A == 0 & data$S == strt[i])
    fitk1 = nnet(data.Y~.,data = datafit[indk1,],size = 2,linout = TRUE,rang = 0.1,decay = 0.02, maxit = 200,maxNwt = 1500)
    fitk0 = nnet(data.Y~.,data = datafit[indk0,],size = 2,linout = TRUE,rang = 0.1,decay = 0.02, maxit = 200,maxNwt = 1500)
    indk = c(indk1,indk0)
    hk1 = predict(fitk1,data$X[indk,])*Y.range + Y.min
    hk0 = predict(fitk0,data$X[indk,])*Y.range + Y.min
    pai[i] = length(indk1)/length(indk)
    R_strat[indk] = R_strat[indk] - ((1-pai[i])*hk1+pai[i]*hk0)
  }
  res = est(data$A,data$S,R,R_strat,strt,strt_num,n)
  return(res)
}

est_neural<-function(data,hidden){
  n = length(data$A)
  strt = unique(data$S)
  strt_num = length(strt)
  pai = numeric(strt_num)
  ind0 = which(data$A == 0)
  ind1 = which(data$A == 1)
  datafit = data.frame(data$Y,data$X)
  Y.range = max(data$Y) - min(data$Y)
  Y.min = min(data$Y)
  datafit$data.Y <- (datafit$data.Y - Y.min)/Y.range
  fit1 = neuralnet(data.Y~.,data = datafit[ind1,],hidden = hidden,linear.output = TRUE)
  fit0 = neuralnet(data.Y~.,data = datafit[ind0,],hidden = hidden,linear.output = TRUE)
  h1 = predict(fit1,data$X)*Y.range + Y.min
  h0 = predict(fit0,data$X)*Y.range + Y.min
  pait = length(which(data$A == 1))/n
  R = data$Y - ((1-pait)*h1+pait*h0)
  R_strat = data$Y
  for(i in 1:strt_num){
    indk1 = which(data$A == 1 & data$S == strt[i])
    indk0 = which(data$A == 0 & data$S == strt[i])
    fitk1 = neuralnet(data.Y~.,data = datafit[indk1,],hidden = hidden,linear.output = TRUE)
    fitk0 = neuralnet(data.Y~.,data = datafit[indk0,],hidden = hidden,linear.output = TRUE)
    indk = c(indk1,indk0)
    hk1 = predict(fitk1,data$X[indk,])*Y.range + Y.min
    hk0 = predict(fitk0,data$X[indk,])*Y.range + Y.min
    pai[i] = length(indk1)/length(indk)
    R_strat[indk] = R_strat[indk] - ((1-pai[i])*hk1+pai[i]*hk0)
  }
  res = est(data$A,data$S,R,R_strat,strt,strt_num,n)
  return(res)
}

est_neural2<-function(data,hidden){
  n = length(data$A)
  strt = unique(data$S)
  strt_num = length(strt)
  pai = numeric(strt_num)
  ind0 = which(data$A == 0)
  ind1 = which(data$A == 1)
  datafit = data.frame(data$Y,data$X)
  fit1 = neuralnet(data.Y~.,data = datafit[ind1,],hidden = hidden,linear.output = TRUE)
  fit0 = neuralnet(data.Y~.,data = datafit[ind0,],hidden = hidden,linear.output = TRUE)
  h1 = predict(fit1,data$X)
  h0 = predict(fit0,data$X)
  pait = length(which(data$A == 1))/n
  R = data$Y - ((1-pait)*h1+pait*h0)
  R_strat = data$Y
  for(i in 1:strt_num){
    indk1 = which(data$A == 1 & data$S == strt[i])
    indk0 = which(data$A == 0 & data$S == strt[i])
    fitk1 = neuralnet(data.Y~.,data = datafit[indk1,],hidden = hidden,linear.output = TRUE)
    fitk0 = neuralnet(data.Y~.,data = datafit[indk0,],hidden = hidden,linear.output = TRUE)
    indk = c(indk1,indk0)
    hk1 = predict(fitk1,data$X[indk,])
    hk0 = predict(fitk0,data$X[indk,])
    pai[i] = length(indk1)/length(indk)
    R_strat[indk] = R_strat[indk] - ((1-pai[i])*hk1+pai[i]*hk0)
  }
  res = est(data$A,data$S,R,R_strat,strt,strt_num,n)
  return(res)
}


est_gbm<-function(data){
  n = length(data$A)
  strt = unique(data$S)
  strt_num = length(strt)
  pai = numeric(strt_num)
  ind0 = which(data$A == 0)
  ind1 = which(data$A == 1)
  datafit = data.frame(data$Y,data$X)
  fit1 = gbm(data.Y~.,data = datafit[ind1,],n.trees = 50,interaction.depth = 3,bag.fraction = 1,distribution = "gaussian")
  fit0 = gbm(data.Y~.,data = datafit[ind0,],n.trees = 50,interaction.depth = 3,bag.fraction = 1,distribution = "gaussian")
  h1 = predict(fit1,data$X)
  h0 = predict(fit0,data$X)
  pait = length(which(data$A == 1))/n
  R = data$Y - ((1-pait)*h1+pait*h0)
  R_strat = data$Y
  for(i in 1:strt_num){
    indk1 = which(data$A == 1 & data$S == strt[i])
    indk0 = which(data$A == 0 & data$S == strt[i])
    fitk1 = gbm(data.Y~.,data = datafit[indk1,],n.trees = 50,interaction.depth = 3,bag.fraction = 1,distribution = "gaussian")
    fitk0 = gbm(data.Y~.,data = datafit[indk0,],n.trees = 50,interaction.depth = 3,bag.fraction = 1,distribution = "gaussian")
    indk = c(indk1,indk0)
    hk1 = predict(fitk1,data$X[indk,])
    hk0 = predict(fitk0,data$X[indk,])
    pai[i] = length(indk1)/length(indk)
    R_strat[indk] = R_strat[indk] - ((1-pai[i])*hk1+pai[i]*hk0)
  }
  res = est(data$A,data$S,R,R_strat,strt,strt_num,n)
  return(res)
}

est_rf<-function(data){
  n = length(data$A)
  strt = unique(data$S)
  strt_num = length(strt)
  pai = numeric(strt_num)
  ind0 = which(data$A == 0)
  ind1 = which(data$A == 1)
  datafit = data.frame(data$Y,data$X)
  fit1 = randomForest(data.Y~.,data = datafit[ind1,],ntree = 500)
  fit0 = randomForest(data.Y~.,data = datafit[ind0,],ntree = 500)
  h1 = predict(fit1,data$X)
  h0 = predict(fit0,data$X)
  pait = length(which(data$A == 1))/n
  R = data$Y - ((1-pait)*h1+pait*h0)
  R_strat = data$Y
  for(i in 1:strt_num){
    indk1 = which(data$A == 1 & data$S == strt[i])
    indk0 = which(data$A == 0 & data$S == strt[i])
    fitk1 = randomForest(data.Y~.,data = datafit[indk1,],ntree=500)
    fitk0 = randomForest(data.Y~.,data = datafit[indk0,],ntree=500)
    indk = c(indk1,indk0)
    hk1 = predict(fitk1,data$X[indk,])
    hk0 = predict(fitk0,data$X[indk,])
    pai[i] = length(indk1)/length(indk)
    R_strat[indk] = R_strat[indk] - ((1-pai[i])*hk1+pai[i]*hk0)
  }
  res = est(data$A,data$S,R,R_strat,strt,strt_num,n)
  return(res)
}


est_glmnet<-function(data){
  n = length(data$A)
  strt = unique(data$S)
  strt_num = length(strt)
  pai = numeric(strt_num)
  ind0 = which(data$A == 0)
  ind1 = which(data$A == 1)
  datafit = data.frame(data$Y,data$X)
  eGrid <- expand.grid(.alpha=seq(0.1,0.3, by=0.1),.lambda=seq(0,1,by=0.01))
  Control <- trainControl(method="cv",p=0.75) 
  fit1 = train(data.Y~.,data = datafit[ind1,],method = "glmnet",tuneGrid=eGrid,trControl = Control)
  fit0 = train(data.Y~.,data = datafit[ind0,],method = "glmnet",tuneGrid=eGrid,trControl = Control)
  h1 = predict(fit1,data$X)
  h0 = predict(fit0,data$X)
  pait = length(which(data$A == 1))/n
  R = data$Y - ((1-pait)*h1+pait*h0)
  R_strat = data$Y
  for(i in 1:strt_num){
    indk1 = which(data$A == 1 & data$S == strt[i])
    indk0 = which(data$A == 0 & data$S == strt[i])
    fitk1 = train(data.Y~.,data = datafit[indk1,],method = "glmnet",tuneGrid=eGrid,trControl = Control)
    fitk0 = train(data.Y~.,data = datafit[indk0,],method = "glmnet",tuneGrid=eGrid,trControl = Control)
    indk = c(indk1,indk0)
    hk1 = predict(fitk1,data$X[indk,])
    hk0 = predict(fitk0,data$X[indk,])
    pai[i] = length(indk1)/length(indk)
    R_strat[indk] = R_strat[indk] - ((1-pai[i])*hk1+pai[i]*hk0)
  }
  res = est(data$A,data$S,R,R_strat,strt,strt_num,n)
  return(res)
}

est_cart<-function(data){
  n = length(data$A)
  strt = unique(data$S)
  strt_num = length(strt)
  pai = numeric(strt_num)
  ind0 = which(data$A == 0)
  ind1 = which(data$A == 1)
  datafit = data.frame(data$Y,data$X)
  tc <- rpart.control(minsplit=20,minbucket=20,maxdepth=4,xval=5,cp=0.005)
  fit1 = rpart(data.Y~.,data = datafit[ind1,],method = "anova",control = tc)
  fit0 = rpart(data.Y~.,data = datafit[ind0,],method = "anova",control = tc)
  h1 = predict(fit1,data$X)
  h0 = predict(fit0,data$X)
  pait = length(which(data$A == 1))/n
  R = data$Y - ((1-pait)*h1+pait*h0)
  R_strat = data$Y
  for(i in 1:strt_num){
    indk1 = which(data$A == 1 & data$S == strt[i])
    indk0 = which(data$A == 0 & data$S == strt[i])
    fitk1 = rpart(data.Y~.,data = datafit[indk1,],method = "anova",control = tc)
    fitk0 = rpart(data.Y~.,data = datafit[indk0,],method = "anova",control = tc)
    indk = c(indk1,indk0)
    hk1 = predict(fitk1,data$X[indk,])
    hk0 = predict(fitk0,data$X[indk,])
    pai[i] = length(indk1)/length(indk)
    R_strat[indk] = R_strat[indk] - ((1-pai[i])*hk1+pai[i]*hk0)
  }
  res = est(data$A,data$S,R,R_strat,strt,strt_num,n)
  return(res)
}

est_cart2<-function(data){
  n = length(data$A)
  strt = unique(data$S)
  strt_num = length(strt)
  pai = numeric(strt_num)
  ind0 = which(data$A == 0)
  ind1 = which(data$A == 1)
  datafit = data.frame(data$Y,data$X)
  tc <- trainControl("cv",10)
  rpart.grid <- expand.grid(.cp=0.005)
  fit1 = train(data.Y~.,data = datafit[ind1,], method="rpart",trControl=tc,tuneGrid=rpart.grid)
  fit0 = train(data.Y~.,data = datafit[ind0,], method="rpart",trControl=tc,tuneGrid=rpart.grid)
  h1 = predict(fit1,data$X)
  h0 = predict(fit0,data$X)
  pait = length(which(data$A == 1))/n
  R = data$Y - ((1-pait)*h1+pait*h0)
  R_strat = data$Y
  for(i in 1:strt_num){
    indk1 = which(data$A == 1 & data$S == strt[i])
    indk0 = which(data$A == 0 & data$S == strt[i])
    fitk1 = train(data.Y~.,data = datafit[indk1,], method="rpart",trControl=tc,tuneGrid=rpart.grid)
    fitk0 = train(data.Y~.,data = datafit[indk0,], method="rpart",trControl=tc,tuneGrid=rpart.grid)
    indk = c(indk1,indk0)
    hk1 = predict(fitk1,data$X[indk,])
    hk0 = predict(fitk0,data$X[indk,])
    pai[i] = length(indk1)/length(indk)
    R_strat[indk] = R_strat[indk] - ((1-pai[i])*hk1+pai[i]*hk0)
  }
  res = est(data$A,data$S,R,R_strat,strt,strt_num,n)
  return(res)
}

ensemble_in<-function(datafit,ind,indpre){
  Y_mean = mean(datafit$data.Y[indpre])
  Y_sd = sd(datafit$data.Y[indpre])
  datafit[indpre,] = as.data.frame(apply(datafit[indpre,],2,scale))

  fit<-caretList(
    data.Y~.,data = datafit[ind,],
    trControl = trainControl(method = "cv",number = 5,savePredictions = "final"),
    methodList = c("lasso","enet","gbm"),
    tuneList = list(
      rf = caretModelSpec(method = "rf",ntree = 200),
      nn = caretModelSpec(method="nnet", tuneGrid = data.frame(.size = 2,.decay = 0.02)),
      cart = caretModelSpec(method = "rpart",tuneGrid = data.frame(.cp=0.005))
    )
  )
  
  fit_en <- caretStack(fit,method="lm",metric="RMSE",
                        trControl=trainControl(method="cv",number=5))
  fit_best <- which.min(c(min(fit[[1]]$results$RMSE),min(fit[[2]]$results$RMSE),min(fit[[3]]$results$RMSE),
                          min(fit[[4]]$results$RMSE),min(fit[[5]]$results$RMSE),min(fit[[6]]$results$RMSE)))
  
  pred_all = predict(fit,newdata = datafit[indpre,])
  pred_en = predict(fit_en,newdata = datafit[indpre,])
  pred_best = pred_all[,fit_best]
  
  return(cbind(pred_all,ensemble = pred_en,best = pred_best)*Y_sd+Y_mean)
}

est_ensemble<-function(data){
  n = length(data$A)
  strt = unique(data$S)
  strt_num = length(strt)
  pai = numeric(strt_num)
  ind0 = which(data$A == 0)
  ind1 = which(data$A == 1)
  pait = length(ind1)/n
  datafit = data.frame(data$Y,data$X)
  
  h1 = ensemble_in(datafit,ind1,1:n)
  h0 = ensemble_in(datafit,ind0,1:n)
  
  n_est <- ncol(h1)
  
  R = matrix(rep(data$Y,n_est),ncol = n_est)
  R_strat = R
  
  R = R - ((1-pait)*h1+pait*h0)
  
  for(i in 1:strt_num){
    indk1 = which(data$A == 1 & data$S == strt[i])
    indk0 = which(data$A == 0 & data$S == strt[i])
    indk = c(indk1,indk0)
    pai[i] = length(indk1)/length(indk)
    
    hk1 = ensemble_in(datafit,indk1,indk)
    hk0 = ensemble_in(datafit,indk0,indk)

    R_strat[indk,] = R_strat[indk,] - ((1-pai[i])*hk1+pai[i]*hk0)
  }
  
  
  res_tau = matrix(0,nrow = 2,ncol = n_est)
  res_var = matrix(0,nrow = 2,ncol = n_est)
  for(j in 1:n_est){
    res_temp = est(data$A,data$S,R[,j],R_strat[,j],strt,strt_num,n)
    res_tau[,j] = res_temp[,1]
    res_var[,j] = res_temp[,2]
  }
  return(list(res_tau = res_tau,res_var = res_var))

}
