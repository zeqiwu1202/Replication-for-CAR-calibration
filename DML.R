samplevar_DML<-function(A,pair,R,R_strat){
  denom = pair*(1-pair)
  tau_est = mean((A-pair)*R/denom)
  tau_strat = mean((A-pair)*R_strat/denom)
  var_est = mean(((A-pair)*R/denom-tau_est)^2)
  var_strat = mean(((A-pair)*R_strat/denom-tau_strat)^2)
  res = matrix(c(tau_est,var_est,tau_strat,var_strat),nrow = 2,ncol = 2,byrow = TRUE)
  return(res)
}

DML_lasso<-function(data){
  n = length(data$A)
  numk = 5
  Folds = createFolds(1:n,k = numk)
  strt = unique(data$S)
  strt_num = length(strt)
  pai = numeric(strt_num)
  ind0 = which(data$A == 0)
  ind1 = which(data$A == 1)
  pait = length(ind1)/n
  data$X = as.matrix(data$X)
  tau_res = matrix(0,nrow = 2, ncol = numk)
  var_res = matrix(0,nrow = 2, ncol = numk)
  R = data$Y
  R_strat = data$Y
  for(k in 1:numk){
    eval(parse(text = paste("indFoldk = Folds$Fold",k,sep = "")))
    indFoldkc = setdiff(1:n,indFoldk)
    fit0<-cv.glmnet(data$X[intersect(ind0,indFoldkc),],data$Y[intersect(ind0,indFoldkc)],family = "gaussian",standardize = FALSE,nfolds = 5)
    fit1<-cv.glmnet(data$X[intersect(ind1,indFoldkc),],data$Y[intersect(ind1,indFoldkc)],family = "gaussian",standardize = FALSE,nfolds = 5)
    beta0<-coef(fit0,s = fit0$lambda.1se)[-1]
    beta1<-coef(fit1,s = fit1$lambda.1se)[-1]
    R[indFoldk] = R[indFoldk] - data$X[indFoldk,]%*%((1-pait)*beta1+pait*beta0)
    for(i in 1:strt_num){
      indk1 = which(data$A == 1 & data$S == strt[i])
      indk0 = which(data$A == 0 & data$S == strt[i])
      fitk0<-cv.glmnet(data$X[intersect(indk0,indFoldkc),],data$Y[intersect(indk0,indFoldkc)],family = "gaussian",standardize = FALSE,nfolds = 5)
      fitk1<-cv.glmnet(data$X[intersect(indk1,indFoldkc),],data$Y[intersect(indk1,indFoldkc)],family = "gaussian",standardize = FALSE,nfolds = 5)
      betak0<-coef(fitk0,s = fitk0$lambda.1se)[-1]
      betak1<-coef(fitk1,s = fitk1$lambda.1se)[-1]
      indk = c(indk1,indk0)
      pai[i] = length(indk1)/length(indk)
      indkpre = intersect(indk,indFoldk)
      R_strat[indkpre] = R_strat[indkpre] - data$X[indkpre,]%*%((1-pai[i])*betak1+pai[i]*betak0)
    }
    pair = sum(data$A[indFoldkc])/length(indFoldkc)
    #res_temp = samplevar_DML(data$A[indFoldk],pair,R[indFoldk],R_strat[indFoldk])
    res_temp = est(data$A[indFoldk],data$S[indFoldk],R[indFoldk],R_strat[indFoldk],
                   strt,strt_num,length(indFoldk))
    tau_res[,k] = res_temp[,1]
    var_res[,k] = res_temp[,2]
  }
  tau_final = apply(tau_res,1,mean)
  var_final = sqrt(apply(var_res^2,1,sum))/5
  return(cbind(tau_final,var_final))
}

DML_rf<-function(data){
  n = length(data$A)
  numk = 5
  Folds = createFolds(1:n,k = numk)
  strt = unique(data$S)
  strt_num = length(strt)
  pai = numeric(strt_num)
  ind0 = which(data$A == 0)
  ind1 = which(data$A == 1)
  pait = length(ind1)/n
  datafit = data.frame(data$Y,data$X)
  tau_res = matrix(0,nrow = 2, ncol = numk)
  var_res = matrix(0,nrow = 2, ncol = numk)
  R = data$Y
  R_strat = data$Y
  for(k in 1:numk){
    eval(parse(text = paste("indFoldk = Folds$Fold",k,sep = "")))
    indFoldkc = setdiff(1:n,indFoldk)
    fit1 = randomForest(data.Y~.,data = datafit[intersect(ind1,indFoldkc),],ntree=50)
    fit0 = randomForest(data.Y~.,data = datafit[intersect(ind0,indFoldkc),],ntree=50)
    h1 = predict(fit1,data$X[indFoldk,])
    h0 = predict(fit0,data$X[indFoldk,])
    pait = length(which(data$A == 1))/n
    R[indFoldk] = R[indFoldk] - ((1-pait)*h1+pait*h0)
    for(i in 1:strt_num){
      indk1 = which(data$A == 1 & data$S == strt[i])
      indk0 = which(data$A == 0 & data$S == strt[i])
      indFoldkc1 <- intersect(indk1,indFoldkc)
      indFoldkc0 <- intersect(indk0,indFoldkc)
      fitk1 = randomForest(data.Y~.,data = datafit[intersect(indk1,indFoldkc),],ntree=50)
      fitk0 = randomForest(data.Y~.,data = datafit[intersect(indk0,indFoldkc),],ntree=50)
      indk = c(indk1,indk0)
      indkpre = intersect(indk,indFoldk)
      hk1 = predict(fitk1,data$X[indkpre,])
      hk0 = predict(fitk0,data$X[indkpre,])
      pai[i] = length(indk1)/length(indk)
      R_strat[indkpre] = R_strat[indkpre] - ((1-pai[i])*hk1+pai[i]*hk0)
    }
    pair <- sum(data$A[indFoldkc])/length(indFoldkc)
    #res_temp_DML = samplevar_DML(data$A[indFoldk],pair,R[indFoldk],R_strat[indFoldk])
    res_temp = est(data$A[indFoldk],data$S[indFoldk],R[indFoldk],R_strat[indFoldk],
                   strt,strt_num,length(indFoldk))
    tau_res[,k] = res_temp[,1]
    var_res[,k] = res_temp[,2]
  }
  tau_final = apply(tau_res,1,mean)
  var_final = sqrt(apply(var_res^2,1,sum))/5
  return(cbind(tau_final,var_final))
}

DML_gbm<-function(data){
  n = length(data$A)
  numk = 5
  Folds = createFolds(1:n,k = numk)
  strt = unique(data$S)
  strt_num = length(strt)
  pai = numeric(strt_num)
  ind0 = which(data$A == 0)
  ind1 = which(data$A == 1)
  pait = length(ind1)/n
  datafit = data.frame(data$Y,data$X)
  tau_res = matrix(0,nrow = 2, ncol = numk)
  var_res = matrix(0,nrow = 2, ncol = numk)
  R = data$Y
  R_strat = data$Y
  for(k in 1:numk){
    eval(parse(text = paste("indFoldk = Folds$Fold",k,sep = "")))
    indFoldkc = setdiff(1:n,indFoldk)
    fit1 = gbm(data.Y~.,data = datafit[intersect(ind1,indFoldkc),],n.trees = 50,interaction.depth = 3,bag.fraction = 1,distribution = "gaussian")
    fit0 = gbm(data.Y~.,data = datafit[intersect(ind0,indFoldkc),],n.trees = 50,interaction.depth = 3,bag.fraction = 1,distribution = "gaussian")
    h1 = predict(fit1,data$X[indFoldk,])
    h0 = predict(fit0,data$X[indFoldk,])
    pait = length(which(data$A == 1))/n
    R[indFoldk] = R[indFoldk] - ((1-pait)*h1+pait*h0)
    for(i in 1:strt_num){
      indk1 = which(data$A == 1 & data$S == strt[i])
      indk0 = which(data$A == 0 & data$S == strt[i])
      fitk1 = gbm(data.Y~.,data = datafit[intersect(indk1,indFoldkc),],n.trees = 50,interaction.depth = 3,bag.fraction = 1,distribution = "gaussian")
      fitk0 = gbm(data.Y~.,data = datafit[intersect(indk0,indFoldkc),],n.trees = 50,interaction.depth = 3,bag.fraction = 1,distribution = "gaussian")
      indk = c(indk1,indk0)
      indkpre = intersect(indk,indFoldk)
      hk1 = predict(fitk1,data$X[indkpre,])
      hk0 = predict(fitk0,data$X[indkpre,])
      pai[i] = length(indk1)/length(indk)
      R_strat[indkpre] = R_strat[indkpre] - ((1-pai[i])*hk1+pai[i]*hk0)
    }
    #res_temp = samplevar_DML(data$A[indFoldk],R[indFoldk],R_strat[indFoldk])
    res_temp = est(data$A[indFoldk],data$S[indFoldk],R[indFoldk],R_strat[indFoldk],
                   strt,strt_num,length(indFoldk))
    tau_res[,k] = res_temp[,1]
    var_res[,k] = res_temp[,2]
  }
  tau_final = apply(tau_res,1,mean)
  var_final = sqrt(apply(var_res^2,1,sum))/5
  return(cbind(tau_final,var_final))
}

DML_nnet<-function(data){
  n = length(data$A)
  numk = 5
  Folds = createFolds(1:n,k = numk)
  strt = unique(data$S)
  strt_num = length(strt)
  pai = numeric(strt_num)
  ind0 = which(data$A == 0)
  ind1 = which(data$A == 1)
  pait = length(ind1)/n
  datafit = data.frame(data$Y,data$X)
  Y.range = max(data$Y) - min(data$Y)
  Y.min = min(data$Y)
  datafit$data.Y <- (datafit$data.Y - Y.min)/Y.range
  tau_res = matrix(0,nrow = 2, ncol = numk)
  var_res = matrix(0,nrow = 2, ncol = numk)
  R = data$Y
  R_strat = data$Y
  for(k in 1:numk){
    eval(parse(text = paste("indFoldk = Folds$Fold",k,sep = "")))
    indFoldkc = setdiff(1:n,indFoldk)
    fit1 = nnet(data.Y~.,data = datafit[intersect(ind1,indFoldkc),],size = 2,linout = TRUE,rang = 0.1,decay = 0.02, maxit = 200,maxNwt = 1500)
    fit0 = nnet(data.Y~.,data = datafit[intersect(ind0,indFoldkc),],size = 2,linout = TRUE,rang = 0.1,decay = 0.02, maxit = 200,maxNwt = 1500)
    h1 = predict(fit1,data$X[indFoldk,])*Y.range + Y.min
    h0 = predict(fit0,data$X[indFoldk,])*Y.range + Y.min
    pait = length(which(data$A == 1))/n
    R[indFoldk] = R[indFoldk] - ((1-pait)*h1+pait*h0)
    for(i in 1:strt_num){
      indk1 = which(data$A == 1 & data$S == strt[i])
      indk0 = which(data$A == 0 & data$S == strt[i])
      fitk1 = nnet(data.Y~.,data = datafit[intersect(indk1,indFoldkc),],size = 2,linout = TRUE,rang = 0.1,rang = 0.1,decay = 0.02, maxit = 200,maxNwt = 1500)
      fitk0 = nnet(data.Y~.,data = datafit[intersect(indk0,indFoldkc),],size = 2,linout = TRUE,rang = 0.1,rang = 0.1,decay = 0.02, maxit = 200,maxNwt = 1500)
      indk = c(indk1,indk0)
      indkpre = intersect(indk,indFoldk)
      hk1 = predict(fitk1,data$X[indkpre,])*Y.range + Y.min
      hk0 = predict(fitk0,data$X[indkpre,])*Y.range + Y.min
      pai[i] = length(indk1)/length(indk)
      R_strat[indkpre] = R_strat[indkpre] - ((1-pai[i])*hk1+pai[i]*hk0)
    }
    #res_temp = samplevar_DML(data$A[indFoldk],R[indFoldk],R_strat[indFoldk])
    res_temp = est(data$A[indFoldk],data$S[indFoldk],R[indFoldk],R_strat[indFoldk],
                   strt,strt_num,length(indFoldk))
    tau_res[,k] = res_temp[,1]
    var_res[,k] = res_temp[,2]
  }
  tau_final = apply(tau_res,1,mean)
  var_final = sqrt(apply(var_res^2,1,sum))/5
  return(cbind(tau_final,var_final))
}

DML_neural<-function(data,hidden){
  n = length(data$A)
  numk = 5
  Folds = createFolds(1:n,k = numk)
  strt = unique(data$S)
  strt_num = length(strt)
  pai = numeric(strt_num)
  ind0 = which(data$A == 0)
  ind1 = which(data$A == 1)
  pait = length(ind1)/n
  datafit = data.frame(data$Y,data$X)
  Y.range = max(data$Y) - min(data$Y)
  Y.min = min(data$Y)
  datafit$data.Y <- (datafit$data.Y - Y.min)/Y.range
  tau_res = matrix(0,nrow = 2, ncol = numk)
  var_res = matrix(0,nrow = 2, ncol = numk)
  R = data$Y
  R_strat = data$Y
  for(k in 1:numk){
    eval(parse(text = paste("indFoldk = Folds$Fold",k,sep = "")))
    indFoldkc = setdiff(1:n,indFoldk)
    fit1 = neuralnet(data.Y~.,data = datafit[intersect(ind1,indFoldkc),],hidden = hidden,linear.output = TRUE)
    fit0 = neuralnet(data.Y~.,data = datafit[intersect(ind0,indFoldkc),],hidden = hidden,linear.output = TRUE)
    h1 = predict(fit1,data$X[indFoldk,])*Y.range + Y.min
    h0 = predict(fit0,data$X[indFoldk,])*Y.range + Y.min
    pait = length(which(data$A == 1))/n
    R[indFoldk] = R[indFoldk] - ((1-pait)*h1+pait*h0)
    for(i in 1:strt_num){
      indk1 = which(data$A == 1 & data$S == strt[i])
      indk0 = which(data$A == 0 & data$S == strt[i])
      fitk1 = neuralnet(data.Y~.,data = datafit[intersect(indk1,indFoldkc),],hidden = hidden,linear.output = TRUE)
      fitk0 = neuralnet(data.Y~.,data = datafit[intersect(indk0,indFoldkc),],hidden = hidden,linear.output = TRUE)
      indk = c(indk1,indk0)
      indkpre = intersect(indk,indFoldk)
      hk1 = predict(fitk1,data$X[indkpre,])*Y.range + Y.min
      hk0 = predict(fitk0,data$X[indkpre,])*Y.range + Y.min
      pai[i] = length(indk1)/length(indk)
      R_strat[indkpre] = R_strat[indkpre] - ((1-pai[i])*hk1+pai[i]*hk0)
    }
    #res_temp = samplevar_DML(data$A[indFoldk],R[indFoldk],R_strat[indFoldk])
    res_temp = est(data$A[indFoldk],data$S[indFoldk],R[indFoldk],R_strat[indFoldk],
                   strt,strt_num,length(indFoldk))
    tau_res[,k] = res_temp[,1]
    var_res[,k] = res_temp[,2]
  }
  tau_final = apply(tau_res,1,mean)
  var_final = sqrt(apply(var_res^2,1,sum))/5
  return(cbind(tau_final,var_final))
}

DML_neural_deep<-function(data){
  n = length(data$A)
  numk = 5
  Folds = createFolds(1:n,k = numk)
  strt = unique(data$S)
  strt_num = length(strt)
  pai = numeric(strt_num)
  ind0 = which(data$A == 0)
  ind1 = which(data$A == 1)
  pait = length(ind1)/n
  datafit = data.frame(data$Y,data$X)
  Y.range = max(data$Y) - min(data$Y)
  Y.min = min(data$Y)
  datafit$data.Y <- (datafit$data.Y - Y.min)/Y.range
  tau_res = matrix(0,nrow = 2, ncol = numk)
  var_res = matrix(0,nrow = 2, ncol = numk)
  R = data$Y
  R_strat = data$Y
  for(k in 1:numk){
    eval(parse(text = paste("indFoldk = Folds$Fold",k,sep = "")))
    indFoldkc = setdiff(1:n,indFoldk)
    fit1 = neuralnet(data.Y~.,data = datafit[intersect(ind1,indFoldkc),],hidden = c(5,2),linear.output = TRUE)
    fit0 = neuralnet(data.Y~.,data = datafit[intersect(ind0,indFoldkc),],hidden = c(5,2),linear.output = TRUE)
    h1 = predict(fit1,data$X[indFoldk,])*Y.range + Y.min
    h0 = predict(fit0,data$X[indFoldk,])*Y.range + Y.min
    pait = length(which(data$A == 1))/n
    R[indFoldk] = R[indFoldk] - ((1-pait)*h1+pait*h0)
    for(i in 1:strt_num){
      indk1 = which(data$A == 1 & data$S == strt[i])
      indk0 = which(data$A == 0 & data$S == strt[i])
      fitk1 = neuralnet(data.Y~.,data = datafit[intersect(indk1,indFoldkc),],hidden = c(5,2),linear.output = TRUE)
      fitk0 = neuralnet(data.Y~.,data = datafit[intersect(indk0,indFoldkc),],hidden = c(5,2),linear.output = TRUE)
      indk = c(indk1,indk0)
      indkpre = intersect(indk,indFoldk)
      hk1 = predict(fitk1,data$X[indkpre,])*Y.range + Y.min
      hk0 = predict(fitk0,data$X[indkpre,])*Y.range + Y.min
      pai[i] = length(indk1)/length(indk)
      R_strat[indkpre] = R_strat[indkpre] - ((1-pai[i])*hk1+pai[i]*hk0)
    }
    #res_temp = samplevar_DML(data$A[indFoldk],R[indFoldk],R_strat[indFoldk])
    res_temp = est(data$A[indFoldk],data$S[indFoldk],R[indFoldk],R_strat[indFoldk],
                   strt,strt_num,length(indFoldk))
    tau_res[,k] = res_temp[,1]
    var_res[,k] = res_temp[,2]
  }
  tau_final = apply(tau_res,1,mean)
  var_final = sqrt(apply(var_res^2,1,sum))/5
  return(cbind(tau_final,var_final))
}

DML_cart<-function(data){
  n = length(data$A)
  numk = 5
  Folds = createFolds(1:n,k = numk)
  strt = unique(data$S)
  strt_num = length(strt)
  pai = numeric(strt_num)
  ind0 = which(data$A == 0)
  ind1 = which(data$A == 1)
  pait = length(ind1)/n
  datafit = data.frame(data$Y,data$X)
  tau_res = matrix(0,nrow = 2, ncol = numk)
  var_res = matrix(0,nrow = 2, ncol = numk)
  R = data$Y
  R_strat = data$Y
  tc <- rpart.control(minsplit=20,minbucket=20,maxdepth=3,xval=5,cp=0.005)
  for(k in 1:numk){
    eval(parse(text = paste("indFoldk = Folds$Fold",k,sep = "")))
    indFoldkc = setdiff(1:n,indFoldk)
    fit1 = rpart(data.Y~.,data = datafit[intersect(ind1,indFoldkc),],method = "anova",control = tc)
    fit0 = rpart(data.Y~.,data = datafit[intersect(ind0,indFoldkc),],method = "anova",control = tc)
    h1 = predict(fit1,data$X[indFoldk,])
    h0 = predict(fit0,data$X[indFoldk,])
    pait = length(which(data$A == 1))/n
    R[indFoldk] = R[indFoldk] - ((1-pait)*h1+pait*h0)
    for(i in 1:strt_num){
      indk1 = which(data$A == 1 & data$S == strt[i])
      indk0 = which(data$A == 0 & data$S == strt[i])
      fitk1 = rpart(data.Y~.,data = datafit[intersect(indk1,indFoldkc),],method = "anova",control = tc)
      fitk0 = rpart(data.Y~.,data = datafit[intersect(indk0,indFoldkc),],method = "anova",control = tc)
      indk = c(indk1,indk0)
      indkpre = intersect(indk,indFoldk)
      hk1 = predict(fitk1,data$X[indkpre,])
      hk0 = predict(fitk0,data$X[indkpre,])
      pai[i] = length(indk1)/length(indk)
      R_strat[indkpre] = R_strat[indkpre] - ((1-pai[i])*hk1+pai[i]*hk0)
    }
    #res_temp = samplevar_DML(data$A[indFoldk],R[indFoldk],R_strat[indFoldk])
    res_temp = est(data$A[indFoldk],data$S[indFoldk],R[indFoldk],R_strat[indFoldk],
                   strt,strt_num,length(indFoldk))
    tau_res[,k] = res_temp[,1]
    var_res[,k] = res_temp[,2]
  }
  tau_final = apply(tau_res,1,mean)
  var_final = sqrt(apply(var_res^2,1,sum))/5
  return(cbind(tau_final,var_final))
}

DML_ensemble<-function(data,name_methods){
  n = length(data$A)
  numk = 5
  Folds = createFolds(1:n,k = numk)
  strt = unique(data$S)
  strt_num = length(strt)
  pai = numeric(strt_num)
  ind0 = which(data$A == 0)
  ind1 = which(data$A == 1)
  pait = length(ind1)/n
  datafit = data.frame(data$Y,data$X)
  for(j in 1:8){
    eval(parse(text = paste(name_methods[j],"_tau = matrix(0,nrow = 2, ncol = numk)",sep = '')))
    eval(parse(text = paste(name_methods[j],"_var = matrix(0,nrow = 2, ncol = numk)",sep = '')))
  }
  R = matrix(rep(data$Y,8),ncol = 8)
  R_strat = R
  for(k in 1:numk){
    eval(parse(text = paste("indFoldk = Folds$Fold",k,sep = "")))
    indFoldkc = setdiff(1:n,indFoldk)
    h1 = ensemble_in(datafit,intersect(ind1,indFoldkc),indFoldk)
    h0 = ensemble_in(datafit,intersect(ind0,indFoldkc),indFoldk)
    pait = length(which(data$A == 1))/n
    R[indFoldk,] = R[indFoldk,] - ((1-pait)*h1+pait*h0)
    
    for(i in 1:strt_num){
      indk1 = which(data$A == 1 & data$S == strt[i])
      indk0 = which(data$A == 0 & data$S == strt[i])
      indk = c(indk1,indk0)
      indkpre = intersect(indk,indFoldk)
      pai[i] = length(indk1)/length(indk)
      hk1 = ensemble_in(datafit,intersect(indk1,indFoldkc),indkpre)
      hk0 = ensemble_in(datafit,intersect(indk0,indFoldkc),indkpre)
      R_strat[indkpre,] = R_strat[indkpre,] - ((1-pai[i])*hk1+pai[i]*hk0)
    }
    pair <- sum(data$A[indFoldkc])/length(indFoldkc)
    for(j in 1:8){
      eval(parse(text = paste(name_methods[j],"_temp = samplevar_DML(data$A[indFoldk],pair,R[indFoldk,j],R_strat[indFoldk,j])",sep = '')))
      eval(parse(text = paste(name_methods[j],"_tau[,k] = ",name_methods[j],"_temp[,1]",sep = '')))
      eval(parse(text = paste(name_methods[j],"_var[,k] = ",name_methods[j],"_temp[,2]",sep = '')))
    }

  }
  
  tau_final = matrix(0,nrow = 2,ncol = 8)
  var_final = matrix(0,nrow = 2,ncol = 8)
  for(m in 1:8){
    eval(parse(text = paste("tau_final[,m] = apply(",name_methods[m],"_tau,1,mean)",sep = '')))
    eval(parse(text = paste("var_final[,m] = apply(",name_methods[m],"_var,1,mean)",sep = '')))
  }
  return(list(res_tau = tau_final,res_var = sqrt(var_final/n)))
}
