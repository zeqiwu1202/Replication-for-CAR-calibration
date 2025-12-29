library(ggplot2)
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


data <- read_dta("Data/Data_Uganda/glaseu_four_rounds.dta")
data_b <- read_dta("Data/Data_Uganda/glaseu_baseline.dta")
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

strata_index = rep(0,length(data$hhid))
stata_dummy <- data[, paste0("strata", 1:78)]
for (i in 1:78){
  strata_index[data[, paste0("strata", i)] == 1] = i
}

data_Malawi = data.frame(Y=data$mv1_amountsaved_resp_tot * 0.005928101, X = data$b_amountsaved_resp_tot2, A=data$treated,S=strata_index)
data_Malawi = na.omit(data_Malawi)








cal_estimates<-function(data, data_information,location){
  
  ind1 = which(data_information$A == 1)
  ind0 = which(data_information$A == 0)
  datafit = data.frame(data_information$Y,data_information$X)
  fit1_rf_malawi = randomForest(data_information.Y~.,data = datafit[ind1,],ntree=40)
  fit0_rf_malawi = randomForest(data_information.Y~.,data = datafit[ind0,],ntree=40)
  
  S_count = table(data$S)
  keep_stata = names(S_count[S_count > 6])
  data <- data[(data$S %in% keep_stata),]
  n = length(data$A)
  print(n)
  print(length(unique(data$S)))
  numk = 1
  Folds = createFolds_by_strat(1:n, numk, data$S, data$A)
  strt = unique(data$S)
  strt_num = length(strt)
  pai = numeric(strt_num)
  ind0 = which(data$A == 0)
  ind1 = which(data$A == 1)
  datafit = data.frame(data$Y,data$X)
  tau_res = matrix(0,nrow = 7, ncol = numk)
  sd_res = matrix(0,nrow = 7, ncol = numk)
  n = length(data$Y)
  hX_eachFold_1_strat = data.frame(zero = rep(0, n), X = rep(0, n), rfmalawi = rep(0, n))
  hX_eachFold_0_strat = data.frame(zero = rep(0, n), logX = rep(0, n), rfmalawi = rep(0, n))
  hX_reg_adjusted_0 = data.frame(reg = rep(0, n))
  hX_reg_adjusted_1 = data.frame(reg = rep(0, n))
  pai_each_data = rep(0, length(data$A))
  str_indicator = rep(0, length(data$A))
  R_dml2 = rep(0, length(data$A))
  R_dml = rep(0, length(data$A))
  R_dml_rf = rep(0, length(data$A))
  R_dml_nn = rep(0, length(data$A))
  R_dml_lm = rep(0, length(data$A))
  R_naive = rep(0,length(data$Y))
  R_0 = rep(0, length(data$A))
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
  for(i in 1:strt_num){
    
    indk1 = which(data$A == 1 & data$S == strt[i])
    indk0 = which(data$A == 0 & data$S == strt[i])
    indk = c(indk1,indk0)
    str_indicator[indk] = i
    pi_for_this_str = 0.5
    pai_each_data[indk] = pi_for_this_str

      pi_for_this_str = length(indk1)/length(indk)
      pai_each_data[indk] = pi_for_this_str
      
      
      # rf from malawi dataset
      datafit_malawi = data.frame(data_information.X = data$X[indk])
      hk1_malawi = predict(fit1_rf_malawi,datafit_malawi)
      hk0_malawi = predict(fit0_rf_malawi,datafit_malawi)
      hX_eachFold_1_strat$rfmalawi[indk] = hk1_malawi
      hX_eachFold_0_strat$rfmalawi[indk] = hk0_malawi
      
      # linear regression
      logy = log(data$Y + 1)
      logx = log(data$X + 1)
      log_lm = lm(logy~logx)
      hX_eachFold_1_strat$X[indk] = data$X[indk]
      hX_eachFold_0_strat$logX[indk] = (data$X[indk]+1) ^ (log_lm$coefficients[2])
    }
  
  
  est_temp_cal_X = cal_est(strt_num,
                           str_indicator,
                           data$A,
                           pai_each_data,
                           subset(hX_eachFold_1_strat,select="X"),
                           subset(hX_eachFold_0_strat,select="zero"),
                           data$Y,
                           population_pi = 0.5) 
  est_temp_cal_logX = cal_est(strt_num,
                              str_indicator,
                              data$A,
                              pai_each_data,
                              subset(hX_eachFold_1_strat,select="zero"),
                              subset(hX_eachFold_0_strat,select="logX"),
                              data$Y,
                              population_pi = 0.5) 
  est_temp_cal_XlogX = cal_est(strt_num,
                               str_indicator,
                               data$A,
                               pai_each_data,
                               subset(hX_eachFold_1_strat,select="X"),
                               subset(hX_eachFold_0_strat,select="logX"),
                               data$Y,
                               population_pi = 0.5) 
  
  
  est_temp_cal_X_malawi = cal_est(strt_num,
                                  str_indicator,
                                  data$A,
                                  pai_each_data,
                                  subset(hX_eachFold_1_strat,select=c("X","rfmalawi")),
                                  subset(hX_eachFold_0_strat,select=c("zero","rfmalawi")),
                                  data$Y,
                                  population_pi = 0.5) 
  est_temp_cal_logX_malawi = cal_est(strt_num,
                                     str_indicator,
                                     data$A,
                                     pai_each_data,
                                     subset(hX_eachFold_1_strat,select=c("zero","rfmalawi")),
                                     subset(hX_eachFold_0_strat,select=c("logX","rfmalawi")),
                                     data$Y,
                                     population_pi = 0.5)
  est_temp_cal_X_logX_malawi = cal_est(strt_num,
                                       str_indicator,
                                       data$A,
                                       pai_each_data,
                                       subset(hX_eachFold_1_strat,select=c("X","rfmalawi")),
                                       subset(hX_eachFold_0_strat,select=c("logX","rfmalawi")),
                                       data$Y,
                                       population_pi = 0.5) 
  
  
  
  
  
  
  
  est_temp_sdim = cal_est(strt_num,
                          str_indicator,
                          data$A,
                          pai_each_data,
                          subset(hX_eachFold_1_strat,select=c("X","rfmalawi")),
                          subset(hX_eachFold_0_strat,select=c("logX","rfmalawi")),
                          data$Y,
                          uniform_weights = TRUE) 
  
  
  
  est_temp_dml_lm = est(data$A,data$S,R_dml_lm,R_dml_lm,
                        strt,strt_num,length(data$A)) # DML_naive_lm
  
  k = 1
  tau_res[1,k] = est_temp_sdim[1]
  sd_res[1,k] = est_temp_sdim[2]
  
  tau_res[2,k] = est_temp_cal_X[1]
  sd_res[2,k] = est_temp_cal_X[2]
  
  tau_res[3,k] = est_temp_cal_logX[1]
  sd_res[3,k] = est_temp_cal_logX[2]
  
  tau_res[4,k] = est_temp_cal_XlogX[1]
  sd_res[4,k] = est_temp_cal_XlogX[2]
  
  tau_res[5,k] = est_temp_cal_X_malawi[1]
  sd_res[5,k] = est_temp_cal_X_malawi[2]
  
  tau_res[6,k] = est_temp_cal_logX_malawi[1]
  sd_res[6,k] = est_temp_cal_logX_malawi[2]
  
  tau_res[7,k] = est_temp_cal_X_logX_malawi[1]
  sd_res[7,k] = est_temp_cal_X_logX_malawi[2]
  
  
  
  
  data_plot <- data.frame(
    estimate = tau_res,        
    lower = tau_res - 1.96 * sd_res,          
    upper = tau_res + 1.96 * sd_res,          
    variable = c("A", "B", "C", "D","E","F","G")         
  )
  
  p<-ggplot(data_plot, aes(x = variable, y = estimate)) +
    geom_point(size = 2) + 
    geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.1) +  
    theme_minimal() +  
    theme(panel.grid.major = element_blank(),  
          panel.grid.minor = element_blank(),
          plot.title=element_text(hjust=0.5),
          panel.border = element_rect(color = "black", fill = NA, size = 1),
          axis.text.x = element_text(color = "black",family = "mono",size=10, vjust = 0.2),
          axis.text.y = element_text(color = "black", size=10),
          axis.title.x = element_blank())+
    geom_hline(yintercept = 0, linetype = "dashed", color = "black", size = 0.5)+
    labs(y = "Estimated ATE", title = paste("The ATE estimates in ", location, sep = ""))+
    scale_x_discrete(labels = c("A" = "sdim", 
                                "B" = expression(paste(cal_X)), 
                                "C" = expression(paste(cal_X^{beta})), 
                                "D" = expression(paste(cal_X_X^{beta})),
                                "E" = expression(paste(cal_info_X)),
                                "F" = expression(paste(cal_info_X^{beta})),
                                "G" = expression(paste(cal_info_X_X^{beta}))
                                ))
  
  title = paste("ATE_estimates_in_", location, ".pdf",sep="")
  ggsave(title, plot = p, device = "pdf", width = 8, height = 2.5)
  
  
  
  result = data.frame(est=round(tau_res,3),sd=round(sd_res,3))
  result = t(result)
  return(result)
}

data = data_Uganda
data_information = data_Malawi
result_Uganda = cal_estimates(data, data_information, "Uganda")
print(result_Uganda,quote=FALSE)
data = data_Malawi
data_information =  data_Uganda 
result_Malawi = cal_estimates(data, data_information, "Malawi")
print(result_Malawi,quote=FALSE)
write.csv(rbind(result_Uganda,result_Malawi), file = "empirical.csv", row.names = FALSE)




