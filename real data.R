library(haven)
library(pls)
library(mda)
library(splines)
library(randomForest)
library(caret)
library(dplyr)
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
library(carat)
Rcpp::sourceCpp("vestimator.cpp")
source("estimators.R")
source("cal_estimator.R")
source("calibration.R")
source("createFolds_by_strat.R")
source("Model1.R")
source("Model2.R")
source("Model3.R")
source("Model4.R")
source("tables.R")
source("DML.R")
source("parallel.R")
RNGkind("L'Ecuyer-CMRG")
set.seed(0)
# mc.reset.stream()


data <- read_dta("Data/Data_Uganda/glaseu_four_rounds.dta")

vars <- grep("_amountsaved_resp_(home|bank|sacco|rosca|friend|mobile|shop|leader|farmgroup)$",
             names(data), value = TRUE)

for (v in vars) {
  cutoff <- quantile(data[[v]], 0.99, na.rm = TRUE)
  data[[v]] <- ifelse(data[[v]] > cutoff & !is.na(data[[v]]), cutoff, data[[v]])
}

for (v in c("b", "mv1", "mv2", "e")) {
  data[[paste0(v, "_amountsaved_resp_ocash2")]] <-
    data[[paste0(v, "_amountsaved_resp_shop")]] +
    data[[paste0(v, "_amountsaved_resp_leader")]] +
    data[[paste0(v, "_amountsaved_resp_farmgroup")]]
  
  data[[paste0(v, "_amountsaved_resp_totw")]] <-
    data[[paste0(v, "_amountsaved_resp_home")]] +
    data[[paste0(v, "_amountsaved_resp_bank")]] +
    data[[paste0(v, "_amountsaved_resp_sacco")]] +
    data[[paste0(v, "_amountsaved_resp_rosca")]] +
    data[[paste0(v, "_amountsaved_resp_friend")]] +
    data[[paste0(v, "_amountsaved_resp_mobile")]] +
    data[[paste0(v, "_amountsaved_resp_ocash2")]]
  
  data[[paste0(v, "_amountsaved_resp_tot")]] <- data[[paste0(v, "_amountsaved_resp_totw")]]
}




strata_index = rep(0,length(data$hhid))
stata_dummy <- data[, paste0("strata", 1:41)]
for (i in 1:41){
  strata_index[data[, paste0("strata", i)] == 1] = i
}




data_Uganda = data.frame(Y=data$mv1_amountsaved_resp_tot * 0.000368169, X = data$b_amountsaved_resp_tot2, A=data$treated, S = strata_index)

data_Uganda = na.omit(data_Uganda)





data <- read_dta("Data/Data_Malawi/glasem_four_rounds.dta")

vars <- grep("_amountsaved_resp_(home|bank|sacco|rosca|friend|mobile|shop|leader|farmgroup)$",
             names(data), value = TRUE)

for (v in vars) {
  cutoff <- quantile(data[[v]], 0.99, na.rm = TRUE)
  data[[v]] <- ifelse(data[[v]] > cutoff & !is.na(data[[v]]), cutoff, data[[v]])
}

for (v in c("b", "mv1", "mv2", "e")) {
  data[[paste0(v, "_amountsaved_resp_ocash2")]] <-
    data[[paste0(v, "_amountsaved_resp_shop")]] +
    data[[paste0(v, "_amountsaved_resp_leader")]] +
    data[[paste0(v, "_amountsaved_resp_farmgroup")]]
  
  data[[paste0(v, "_amountsaved_resp_totw")]] <-
    data[[paste0(v, "_amountsaved_resp_home")]] +
    data[[paste0(v, "_amountsaved_resp_bank")]] +
    data[[paste0(v, "_amountsaved_resp_sacco")]] +
    data[[paste0(v, "_amountsaved_resp_rosca")]] +
    data[[paste0(v, "_amountsaved_resp_friend")]] +
    data[[paste0(v, "_amountsaved_resp_mobile")]] +
    data[[paste0(v, "_amountsaved_resp_ocash2")]]
  
  data[[paste0(v, "_amountsaved_resp_tot")]] <- data[[paste0(v, "_amountsaved_resp_totw")]]
}

data_Malawi = data.frame(Y=data$mv1_amountsaved_resp_tot * 0.005928101, X = data$b_amountsaved_resp_tot2, A=data$treated)
data_Malawi = na.omit(data_Malawi)
ind1 = which(data_Malawi$A == 1)
ind0 = which(data_Malawi$A == 0)
datafit = data.frame(data_Malawi$Y,data_Malawi$X)
fit1_rf_malawi = randomForest(data_Malawi.Y~.,data = datafit[ind1,],ntree=50)
fit0_rf_malawi = randomForest(data_Malawi.Y~.,data = datafit[ind0,],ntree=50)








data = data_Uganda
n = length(data$A)
numk = 2
Folds = createFolds_by_strat(1:n, numk, data$S, data$A)
strt = unique(data$S)
strt_num = length(strt)
pai = numeric(strt_num)
ind0 = which(data$A == 0)
ind1 = which(data$A == 1)
datafit = data.frame(data$Y,data$X)
tau_res = matrix(0,nrow = 11, ncol = numk)
sd_res = matrix(0,nrow = 11, ncol = numk)
n = length(data$Y)
hX_eachFold_1_strat = data.frame(rf = rep(0, n), nn = rep(0, n), lm = rep(0, n), rfmalawi = rep(0, n))
hX_eachFold_0_strat = data.frame(rf = rep(0, n), nn = rep(0, n), lm = rep(0, n), rfmalawi = rep(0, n))
hX_reg_adjusted_0 = data.frame(reg = rep(0, n))
hX_reg_adjusted_1 = data.frame(reg = rep(0, n))
#hX_eachFold_1_nostrat = data.frame(rf = rep(NA, n))
#hX_eachFold_0_nostrat = data.frame(rf = rep(NA, n))
pai_each_data = rep(0, length(data$A))
str_indicator = rep(0, length(data$A))
R_dml2 = rep(0, length(data$A))
R_dml = rep(0, length(data$A))
R_dml_rf = rep(0, length(data$A))
R_dml_nn = rep(0, length(data$A))
R_dml_lm = rep(0, length(data$A))
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
    str_indicator[indkpre] = i
    pi_for_this_str = 0.5
    pai_each_data[indkpre] = pi_for_this_str
    if (length(indFoldkc1) > 30 & length(indFoldkc0) > 30){
      pi_for_this_str = length(indk1)/length(indk)
      pai_each_data[indkpre] = pi_for_this_str
    
    
    # rf from malawi dataset
    datafit_malawi = data.frame(data_Malawi.X = data$X[indkpre])
    hk1_malawi = predict(fit1_rf_malawi,datafit_malawi)
    hk0_malawi = predict(fit0_rf_malawi,datafit_malawi)
    hX_eachFold_1_strat$rfmalawi[indkpre] = hk1_malawi
    hX_eachFold_0_strat$rfmalawi[indkpre] = hk0_malawi
    
    # linear regression
    
    fitk1 = lm(data.Y~.,data = datafit[intersect(indk1,indFoldkc),])
    fitk0 = lm(data.Y~.,data = datafit[intersect(indk0,indFoldkc),])
    hk1_lm = predict(fitk1,datafit[indkpre,])
    hk0_lm = predict(fitk0,datafit[indkpre,])
    hX_eachFold_1_strat$lm[indkpre] = hk1_lm
    hX_eachFold_0_strat$lm[indkpre] = hk0_lm

    
    
    fitk1 = randomForest(data.Y~.,data = datafit[intersect(indk1,indFoldkc),],ntree=50)
    fitk0 = randomForest(data.Y~.,data = datafit[intersect(indk0,indFoldkc),],ntree=50)
    hk1_rf = predict(fitk1,datafit[indkpre,])
    hk0_rf = predict(fitk0,datafit[indkpre,])
    
    hX_eachFold_1_strat$rf[indkpre] = hk1_rf
    hX_eachFold_0_strat$rf[indkpre] = hk0_rf
    
    
    
    hk1_rf_group = predict(fitk1,datafit[indFoldk,])
    hk0_rf_group = predict(fitk0,datafit[indFoldk,])
    
    hX_eachFold_1_group[[paste("rf", i, sep = "_")]][indFoldk] = hk1_rf_group
    hX_eachFold_0_group[[paste("rf", i, sep = "_")]][indFoldk] = hk0_rf_group
    
    
    fitk1_nn = neuralnet(data.Y~.,data = datafit_nn[intersect(indk1,indFoldkc),],
                         hidden = 2,linear.output = TRUE)
    fitk0_nn = neuralnet(data.Y~.,data = datafit_nn[intersect(indk0,indFoldkc),],
                         hidden = 2,linear.output = TRUE)
    
    
    hk1_nn = predict(fitk1_nn,datafit_nn[indkpre,])*Y.range + Y.min
    hk0_nn = predict(fitk0_nn,datafit_nn[indkpre,])*Y.range + Y.min
    
    hX_eachFold_1_strat$nn[indkpre] = hk1_nn
    hX_eachFold_0_strat$nn[indkpre] = hk0_nn
    
    
    
    hk1_nn_group = predict(fitk1_nn,datafit_nn[indFoldk,])*Y.range + Y.min
    hk0_nn_group = predict(fitk0_nn,datafit_nn[indFoldk,])*Y.range + Y.min
    
    hX_eachFold_1_group[[paste("nn", i, sep = "_")]][indFoldk] = hk1_nn_group
    hX_eachFold_0_group[[paste("nn", i, sep = "_")]][indFoldk] = hk0_nn_group
    
    
    
    R_dml_rf[indkpre] = R_dml[indkpre] - ((1-pi_for_this_str)*hk1_rf+pi_for_this_str*hk0_rf)
    R_dml_nn[indkpre] = R_dml[indkpre] - ((1-pi_for_this_str)*hk1_nn+pi_for_this_str*hk0_nn)
    R_dml_lm[indkpre] = R_dml[indkpre] - ((1-pi_for_this_str)*hk1_lm+pi_for_this_str*hk0_lm)
    }
  }
  # print(max(is.na(hX_eachFold_0_strat[indFoldk,])))
  # print(max(is.na(hX_eachFold_1_strat[indFoldk,])))
  # print(hX_eachFold_1_strat)
  
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
  est_temp_cal_lm_malawi = cal_est(strt_num,
                                   str_indicator[indFoldk],
                                   data$A[indFoldk],
                                   pai_each_data[indFoldk],
                                   subset(hX_eachFold_1_strat[indFoldk,],select=c("lm","rfmalawi")),
                                   subset(hX_eachFold_0_strat[indFoldk,],select=c("lm","rfmalawi")),
                                   data$Y[indFoldk]) 
  est_temp_cal_rf_malawi = cal_est(strt_num,
                                   str_indicator[indFoldk],
                                   data$A[indFoldk],
                                   pai_each_data[indFoldk],
                                   subset(hX_eachFold_1_strat[indFoldk,],select=c("rf","rfmalawi")),
                                   subset(hX_eachFold_0_strat[indFoldk,],select=c("rf","rfmalawi")),
                                   data$Y[indFoldk]) 
  
  est_temp_cal_rflm_malawi = cal_est(strt_num,
                              str_indicator[indFoldk],
                              data$A[indFoldk],
                              pai_each_data[indFoldk],
                              subset(hX_eachFold_1_strat[indFoldk,],select=c("rf","lm","rfmalawi")),
                              subset(hX_eachFold_0_strat[indFoldk,],select=c("rf","lm","rfmalawi")),
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
  
  
  tau_res[1,k] = est_temp_cal_rf[1]
  sd_res[1,k] = est_temp_cal_rf[2]
  tau_res[2,k] = est_temp_cal_nn[1]
  sd_res[2,k] = est_temp_cal_nn[2]
  tau_res[3,k] = est_temp_cal_rfnn[1]
  sd_res[3,k] = est_temp_cal_rfnn[2]
  
  tau_res[4,k] = est_temp_cal_rflm[1]
  sd_res[4,k] = est_temp_cal_rflm[2]
  
  tau_res[5,k] = est_temp_cal_lm_malawi[1]
  sd_res[5,k] = est_temp_cal_lm_malawi[2]
  tau_res[6,k] = est_temp_cal_rf_malawi[1]
  sd_res[6,k] = est_temp_cal_rf_malawi[2]
  tau_res[7,k] = est_temp_cal_rflm_malawi[1]
  sd_res[7,k] = est_temp_cal_rflm_malawi[2]
  
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

