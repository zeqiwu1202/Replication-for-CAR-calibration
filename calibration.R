source("cal_estimator.R")

sim_fun<-function(Model,method){
  return(match.fun(paste(Model,"_",method,sep = '')))
}

calculate_CP<-function(true.value, estimated.values, estimated.sd){
  critical_p = qnorm(0.975)
  lower.bound = estimated.values - critical_p * estimated.sd
  upper.bound = estimated.values + critical_p * estimated.sd
  coverage <- (lower.bound <= true.value & upper.bound >= true.value)
  return(mean(coverage))
}





lin_adj<-function(data){
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
  res = est(data$A,data$S,R_strat,R_strat,strt,strt_num,n)
  return(res[1,])
}

cal_rflm<-function(data){
  n = length(data$A)
  strt = unique(data$S)
  strt_num = length(strt)
  str_indicator = data$Y
  pai_each_data = data$Y
  hX_eachFold_1_strat = data.frame(rf = rep(NA, n),lm = rep(NA, n))
  hX_eachFold_0_strat = data.frame(rf = rep(NA, n),lm = rep(NA, n))
  datafit = data.frame(data$Y,data$Xbeta)
  ind0 = which(data$A == 0)
  ind1 = which(data$A == 1)
  for(i in 1:strt_num){
    indk1 = which(data$A == 1 & data$S == strt[i])
    indk0 = which(data$A == 0 & data$S == strt[i])
    fitk1 = lm(data.Y~.,data = datafit[indk1,])
    fitk0 = lm(data.Y~.,data = datafit[indk0,])
    indk = c(indk1,indk0)
    hk1 = predict(fitk1,datafit[indk,])
    hk0 = predict(fitk0,datafit[indk,])
    hX_eachFold_1_strat$lm[indk] = hk1
    hX_eachFold_0_strat$lm[indk] = hk0

    
    fitk1 = randomForest(data.Y~.,data = datafit[indk1,],ntree=50)
    fitk0 = randomForest(data.Y~.,data = datafit[indk0,],ntree=50)
    hk1_rf = predict(fitk1,datafit[indk,])
    hk0_rf = predict(fitk0,datafit[indk,])
    hX_eachFold_1_strat$rf[indk] = hk1_rf
    hX_eachFold_0_strat$rf[indk] = hk0_rf
    
    pi_for_this_str = length(indk1)/length(indk)
    pai_each_data[indk] = pi_for_this_str 
    str_indicator[indk] = i
  }
  est_temp_cal_rflm = cal_est(strt_num,
                            str_indicator,
                            data$A,
                            pai_each_data,
                            hX_eachFold_0_strat,
                            hX_eachFold_1_strat,
                            data$Y) 
  return(est_temp_cal_rflm)
}

sdim<-function(data){
  n = length(data$A)
  strt = unique(data$S)
  strt_num = length(strt)
  str_indicator = data$Y
  pai_each_data = data$Y
  ind0 = which(data$A == 0)
  ind1 = which(data$A == 1)
  for(i in 1:strt_num){
    indk1 = which(data$A == 1 & data$S == strt[i])
    indk0 = which(data$A == 0 & data$S == strt[i])
    indk = c(indk1,indk0)
    pi_for_this_str = length(indk1)/length(indk)
    pai_each_data[indk] = pi_for_this_str 
    str_indicator[indk] = i
  }
  est_sdim = cal_est(strt_num,
                     str_indicator,
                     data$A,
                     pai_each_data,
                     data$Y,
                     data$Y,
                     data$Y,
                     uniform_weights = TRUE) 
  return(est_sdim)
}












cal_main<-function(data){
  n = length(data$A)
  numk = 2
  Folds = createFolds_by_strat(1:n, numk, data$S, data$A)
  strt = unique(data$S)
  strt_num = length(strt)
  pai = numeric(strt_num)
  ind0 = which(data$A == 0)
  ind1 = which(data$A == 1)
  pait = length(ind1)/n
  datafit = data.frame(data$Y,data$X)
  tau_res = matrix(0,nrow = 11, ncol = numk)
  sd_res = matrix(0,nrow = 11, ncol = numk)
  n = length(data$Y)
  hX_eachFold_1_strat = data.frame(rf = rep(NA, n), nn = rep(NA, n), lm = rep(NA, n), qr=rep(NA,n))
  hX_eachFold_0_strat = data.frame(rf = rep(NA, n), nn = rep(NA, n), lm = rep(NA, n), qr=rep(NA,n))
  hX_reg_adjusted_0 = data.frame(reg = rep(NA, n))
  hX_reg_adjusted_1 = data.frame(reg = rep(NA, n))
  #hX_eachFold_1_nostrat = data.frame(rf = rep(NA, n))
  #hX_eachFold_0_nostrat = data.frame(rf = rep(NA, n))
  pai_each_data = data$Y
  str_indicator = data$A
  R_dml2 = data$Y
  R_dml = data$Y 
  R_dml_rf = data$Y
  R_dml_nn = data$Y
  R_dml_lm = data$Y
  R_naive = rep(0,length(data$Y))
  # for NN
  datafit_nn = data.frame(data$Y,data$X)
  Y.range = max(data$Y) - min(data$Y)
  Y.min = min(data$Y)
  datafit_nn$data.Y <- (datafit_nn$data.Y - Y.min)/Y.range
  
  colnames_rf <- paste("rf", 1:strt_num, sep = "_")
  colnames_nn <- paste("nn", 1:strt_num, sep = "_")  
  # colnames_lm <- paste("lm", 1:strt_num, sep = "_")  
  hX_eachFold_1_group <- data.frame(matrix(ncol = 2 * strt_num, nrow = n))
  hX_eachFold_0_group <- data.frame(matrix(ncol = 2 * strt_num, nrow = n))
  colnames_all <- c(colnames_rf,colnames_nn) 
  colnames(hX_eachFold_1_group) <- colnames_all
  colnames(hX_eachFold_0_group) <- colnames_all
  for(k in 1:numk){
    eval(parse(text = paste("indFoldk = Folds$Fold",k,sep = "")))
    indFoldkc = setdiff(1:n,indFoldk)
    # fit1 = randomForest(data.Y~.,data = datafit[intersect(ind1,indFoldkc),],ntree=50)
    # fit0 = randomForest(data.Y~.,data = datafit[intersect(ind0,indFoldkc),],ntree=50)
    # hX_eachFold_1_nostrat[indFoldk,1] = predict(fit1,data$X[indFoldk,])
    # hX_eachFold_0_nostrat[indFoldk,1] = predict(fit0,data$X[indFoldk,])
    for(i in 1:strt_num){
      indk1 = which(data$A == 1 & data$S == strt[i])
      indk0 = which(data$A == 0 & data$S == strt[i])
      indFoldkc1 <- intersect(indk1,indFoldkc)
      indFoldkc0 <- intersect(indk0,indFoldkc)
      indk = c(indk1,indk0)
      indkpre = intersect(indk,indFoldk)
      
      # linear regression

      fitk1 = lm(data.Y~.,data = datafit[intersect(indk1,indFoldkc),])
      fitk0 = lm(data.Y~.,data = datafit[intersect(indk0,indFoldkc),])
      hk1_lm = predict(fitk1,datafit[indkpre,])
      hk0_lm = predict(fitk0,datafit[indkpre,])
      hX_eachFold_1_strat$lm[indkpre] = hk1_lm
      hX_eachFold_0_strat$lm[indkpre] = hk0_lm
      
      # quantile regression
      # fitk1 = rq(data.Y~.,data = datafit[intersect(indk1,indFoldkc),],tau=0.5)
      # fitk0 = rq(data.Y~.,data = datafit[intersect(indk0,indFoldkc),],tau=0.5)
      # hk1_qr = predict(fitk1,datafit[indkpre,])
      # hk0_qr = predict(fitk0,datafit[indkpre,])
      # hX_eachFold_1_strat$qr[indkpre] = hk1_qr
      # hX_eachFold_0_strat$qr[indkpre] = hk0_qr
      # 
      
      
      # random forest
      fitk1 = randomForest(data.Y~.,data = datafit[intersect(indk1,indFoldkc),],ntree=50)
      fitk0 = randomForest(data.Y~.,data = datafit[intersect(indk0,indFoldkc),],ntree=50)
      hk1_rf = predict(fitk1,data$X[indkpre,])
      hk0_rf = predict(fitk0,data$X[indkpre,])
      
      hX_eachFold_1_strat$rf[indkpre] = hk1_rf
      hX_eachFold_0_strat$rf[indkpre] = hk0_rf
      
      hk1_rf_group = predict(fitk1,data$X[indFoldk,])
      hk0_rf_group = predict(fitk0,data$X[indFoldk,])
      
      hX_eachFold_1_group[[paste("rf", i, sep = "_")]][indFoldk] = hk1_rf_group
      hX_eachFold_0_group[[paste("rf", i, sep = "_")]][indFoldk] = hk0_rf_group
      
      # nn
      
      fitk1_nn = neuralnet(data.Y~.,data = datafit_nn[intersect(indk1,indFoldkc),],
                           hidden = 2,linear.output = TRUE)
      fitk0_nn = neuralnet(data.Y~.,data = datafit_nn[intersect(indk0,indFoldkc),],
                           hidden = 2,linear.output = TRUE)

      
      hk1_nn = predict(fitk1_nn,data$X[indkpre,])*Y.range + Y.min
      hk0_nn = predict(fitk0_nn,data$X[indkpre,])*Y.range + Y.min
      
      hX_eachFold_1_strat$nn[indkpre] = hk1_nn
      hX_eachFold_0_strat$nn[indkpre] = hk0_nn
      
      
      
      hk1_nn_group = predict(fitk1_nn,data$X[indFoldk,])*Y.range + Y.min
      hk0_nn_group = predict(fitk0_nn,data$X[indFoldk,])*Y.range + Y.min
      
      hX_eachFold_1_group[[paste("nn", i, sep = "_")]][indFoldk] = hk1_nn_group
      hX_eachFold_0_group[[paste("nn", i, sep = "_")]][indFoldk] = hk0_nn_group
      
      
      
      str_indicator[indkpre] = i
      # here we have used the whole dataset
      pi_for_this_str = length(indk1)/length(indk)
      pai_each_data[indkpre] = pi_for_this_str
      
      # data.fit.reg1 = data.frame(Y=data$Y,hX_eachFold_1_strat)
      # data.fit.reg0 = data.frame(Y=data$Y,hX_eachFold_0_strat)
      # reg_fit1 = lm(Y~.-1, data=data.fit.reg1[intersect(indk1,indFoldk),])
      # reg_fit0 = lm(Y~.-1, data=data.fit.reg0[intersect(indk0,indFoldk),])
      # hk1_reg = predict(reg_fit1,hX_eachFold_1_strat[indkpre,,drop=FALSE])
      # hk0_reg = predict(reg_fit0,hX_eachFold_0_strat[indkpre,,drop=FALSE])
      # hX_reg_adjusted_0$reg[indkpre] = hk0_reg
      # hX_reg_adjusted_1$reg[indkpre] = hk1_reg
      # estimate pai from each fold or from the whole dataset?
      # R_dml2[indkpre] = R_dml2[indkpre] - ((1-pi_for_this_str)*hk1_reg+pi_for_this_str*hk0_reg)
      
      
      R_dml_rf[indkpre] = R_dml[indkpre] - ((1-pi_for_this_str)*hk1_rf+pi_for_this_str*hk0_rf)
      R_dml_nn[indkpre] = R_dml[indkpre] - ((1-pi_for_this_str)*hk1_nn+pi_for_this_str*hk0_nn)
      R_dml_lm[indkpre] = R_dml[indkpre] - ((1-pi_for_this_str)*hk1_lm+pi_for_this_str*hk0_lm)
    }
    
    est_temp_cal_rf = cal_est(strt_num,
                             str_indicator[indFoldk],
                             data$A[indFoldk],
                             pai_each_data[indFoldk],
                             subset(hX_eachFold_1_strat[indFoldk,],select="rf"),
                             subset(hX_eachFold_0_strat[indFoldk,],select="rf"),
                             data$Y[indFoldk]) 
    est_temp_cal_nn = cal_est(strt_num,
                              str_indicator[indFoldk],
                              data$A[indFoldk],
                              pai_each_data[indFoldk],
                              subset(hX_eachFold_1_strat[indFoldk,],select="nn"),
                              subset(hX_eachFold_0_strat[indFoldk,],select="nn"),
                              data$Y[indFoldk])
    est_temp_cal_rfnn = cal_est(strt_num,
                              str_indicator[indFoldk],
                              data$A[indFoldk],
                              pai_each_data[indFoldk],
                              subset(hX_eachFold_1_strat[indFoldk,],select=c("rf","nn")),
                              subset(hX_eachFold_0_strat[indFoldk,],select=c("rf","nn")),
                              data$Y[indFoldk])
    est_temp_cal_rflm = cal_est(strt_num,
                                str_indicator[indFoldk],
                                data$A[indFoldk],
                                pai_each_data[indFoldk],
                                subset(hX_eachFold_1_strat[indFoldk,],select=c("rf","lm")),
                                subset(hX_eachFold_0_strat[indFoldk,],select=c("rf","lm")),
                                data$Y[indFoldk])
    est_temp_sdim = cal_est(strt_num,
                            str_indicator[indFoldk],
                            data$A[indFoldk],
                            pai_each_data[indFoldk],
                            subset(hX_eachFold_1_strat[indFoldk,],select=c("rf","lm")),
                            subset(hX_eachFold_0_strat[indFoldk,],select=c("rf","lm")),
                            data$Y[indFoldk],
                            uniform_weights = TRUE)
    
    est_temp_dml_rf = est(data$A[indFoldk],data$S[indFoldk],R_dml_rf[indFoldk],R_dml_rf[indFoldk],
                       strt,strt_num,length(indFoldk)) # DML_naive_rf
    est_temp_dml_nn = est(data$A[indFoldk],data$S[indFoldk],R_dml_nn[indFoldk],R_dml_nn[indFoldk],
                          strt,strt_num,length(indFoldk)) # DML_naive_nn
    est_temp_dml_lm = est(data$A[indFoldk],data$S[indFoldk],R_dml_lm[indFoldk],R_dml_lm[indFoldk],
                          strt,strt_num,length(indFoldk)) # DML_naive_lm
    
    est_temp_cal_rf_group = cal_est(strt_num,
                                    str_indicator[indFoldk],
                                    data$A[indFoldk],
                                    pai_each_data[indFoldk],
                                    hX_eachFold_1_group[indFoldk, grep("rf", colnames_all)],
                                    hX_eachFold_0_group[indFoldk, grep("rf", colnames_all)],
                                    data$Y[indFoldk]) #cal_rf_group
    est_temp_cal_nn_group = cal_est(strt_num,
                                    str_indicator[indFoldk],
                                    data$A[indFoldk],
                                    pai_each_data[indFoldk],
                                    hX_eachFold_1_group[indFoldk, grep("nn", colnames_all)],
                                    hX_eachFold_0_group[indFoldk, grep("nn", colnames_all)],
                                    data$Y[indFoldk]) #cal_nn_group
    est_temp_cal_lm_EL = cal_est(strt_num,
                                    str_indicator[indFoldk],
                                    data$A[indFoldk],
                                    pai_each_data[indFoldk],
                                    subset(hX_eachFold_1_strat[indFoldk,],select=c("lm")),
                                    subset(hX_eachFold_0_strat[indFoldk,],select=c("lm")),
                                    data$Y[indFoldk],EL=TRUE) #cal_lm_qr
    
    
    
    
    
    
    
    tau_res[1,k] = est_temp_cal_rf[1]
    sd_res[1,k] = est_temp_cal_rf[2]
    tau_res[2,k] = est_temp_cal_nn[1]
    sd_res[2,k] = est_temp_cal_nn[2]
    tau_res[3,k] = est_temp_cal_rfnn[1]
    sd_res[3,k] = est_temp_cal_rfnn[2]
    
    tau_res[4,k] = est_temp_cal_rflm[1]
    sd_res[4,k] = est_temp_cal_rflm[2]
    
    tau_res[5,k] = est_temp_cal_rf_group[1]
    sd_res[5,k] = est_temp_cal_rf_group[2]
    tau_res[6,k] = est_temp_cal_nn_group[1]
    sd_res[6,k] = est_temp_cal_nn_group[2]
    tau_res[7,k] = est_temp_cal_lm_EL[1]
    sd_res[7,k] = est_temp_cal_lm_EL[2]
    
    tau_res[8,k] = est_temp_dml_rf[1,1]
    sd_res[8,k] = est_temp_dml_rf[1,2]
    
    tau_res[9,k] = est_temp_dml_nn[1,1]
    sd_res[9,k] = est_temp_dml_nn[1,2]
    
    tau_res[10,k] = est_temp_dml_lm[1,1]
    sd_res[10,k] = est_temp_dml_lm[1,2]
    
    tau_res[11,k] = est_temp_sdim[1]
    sd_res[11,k] = est_temp_sdim[2]
    
  }
  sd_final = sqrt(apply(sd_res^2,1,mean)) / sqrt(numk)
  tau_final = apply(tau_res,1,mean)
  return(data.frame(est=tau_final,sd=sd_final))
}






