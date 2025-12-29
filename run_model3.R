library(quantreg)
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
RNGkind("L'Ecuyer-CMRG")
set.seed(0)
# mc.reset.stream()

run_simulation<- function(Model, n, randomization_method){

alpha0 = 5
alpha1 = 2
betavec0 = c(42,83)
betavec1 = c(30,75)
FUN = sim_fun(Model,randomization_method)
lambda = 0.75
weight = 1
sigma1 = 3
sigma0 = 1
pi = 0.5
p=30
S=300



sd_cal_rf = rep(0,S)
sd_cal_nn = rep(0,S)
sd_cal_rfnn = rep(0,S)
sd_DML_nn = rep(0,S)
sd_DML_rf = rep(0,S)
sd_sdim = rep(0,S)


est_cal_rf = rep(0,S)
est_cal_nn = rep(0,S)
est_cal_rfnn = rep(0,S)
est_DML_nn = rep(0,S)
est_DML_rf = rep(0,S)
est_sdim = rep(0,S)


est_dim = rep(0,S)
sd_dim = rep(0,S)

est_cal_rf_g = rep(0,S)
sd_cal_rf_g = rep(0,S)

est_cal_nn_g = rep(0,S)
sd_cal_nn_g = rep(0,S)

est_cal_lm_EL = rep(0,S)
sd_cal_lm_EL = rep(0,S)

est_lin = rep(0,S)
sd_lin = rep(0,S)

est_cal_rflm = rep(0,S)
sd_cal_rflm = rep(0,S)

for (s in 1:S){
  print(s)
  if (randomization_method == "PS") {
    data = FUN(n,alpha0,alpha1,betavec0,betavec1,pi,sigma0,sigma1,p, weight,lambda)
  }
  else {
    data = FUN(n,alpha0,alpha1,betavec0,betavec1,pi,sigma0,sigma1,p)
  }
  
  
  result = cal_main(data)
  sd_cal_rf[s] = result$sd[1]
  sd_cal_nn[s] = result$sd[2]
  sd_cal_rfnn[s] = result$sd[3]
  est_cal_rf[s] = result$est[1]
  est_cal_nn[s] = result$est[2]
  est_cal_rfnn[s] = result$est[3]
  
  

  
  est_cal_rflm[s] = result$est[4]
  sd_cal_rflm[s] = result$sd[4]


  
  
  
  est_cal_rf_g[s] = result$est[5]
  sd_cal_rf_g[s] = result$sd[5]
  
  est_cal_nn_g[s] = result$est[6]
  sd_cal_nn_g[s] = result$sd[6]
  
  est_cal_lm_EL[s] = result$est[7]
  sd_cal_lm_EL[s] = result$sd[7]
  
  est_DML_rf[s] = result$est[8]
  sd_DML_rf[s] = result$sd[8]
  
  est_DML_nn[s] = result$est[9]
  sd_DML_nn[s] = result$sd[9]
 
  
  
  est_lin[s] = result$est[10]
  sd_lin[s] = result$sd[10]
  
  
  est_sdim[s] = result$est[11]
  sd_sdim[s] = result$sd[11]
  
  result_sdim = sdim(data)
  est_sdim[s] = result_sdim[1]
  sd_sdim[s] = result_sdim[2]
  
  # result_rflm = cal_rflm(data)
  # est_cal_rflm[s] = result_rflm[1]
  # sd_cal_rflm[s] = result_rflm[2]
  # result_lin = lin_adj(data)
  # est_lin[s] = result_lin[1]
  # sd_lin[s] = result_lin[2]
  
}


true_value = trueval(Model,alpha0,alpha1,betavec0,betavec1)
data_list <- list(est_cal_rf, est_cal_nn,est_cal_rfnn,est_cal_rflm,
                  est_cal_rf_g,est_cal_nn_g,est_cal_lm_EL,
                  est_DML_rf,est_DML_nn, 
                  est_lin,est_sdim)
sd_data_list <- list(sd_cal_rf, sd_cal_nn,sd_cal_rfnn,sd_cal_rflm,
                     sd_cal_rf_g,sd_cal_nn_g,sd_cal_lm_EL,
                     sd_DML_rf, sd_DML_nn,
                     sd_lin,sd_sdim)
names = c("cal_rf", "cal_nn","cal_rfnn","cal_rflin","cal_rf_g","cal_nn_g","cal_lm_EL",
          "AIPW_rf" , "AIPW_nn","lin","sdim")

Q1 <- lapply(data_list,  function(x) quantile(x, 0.25))
Q3 <- lapply(data_list,  function(x) quantile(x, 0.75))
Q1 = unlist(Q1)
Q3 = unlist(Q3)
IQR <- Q3 - Q1

lower_bound <- min(Q1 - 2 * IQR)
upper_bound <- max(Q3 + 2 * IQR)

bounds = c(lower_bound,upper_bound)
file_name = sprintf("20251008_%s_n%s_%s.pdf", Model, n, randomization_method)
pdf(file_name, width = 12, height = 8)

boxplot(data_list, names=names, main = "Comparison of the Methods",
        col = c("skyblue", "lightgreen", "lightpink", "lightyellow", "lightcoral","lightsalmon","lavender",
                "lightseagreen","lightsteelblue","peachpuff","thistle"),
        ylim=bounds)

abline(h = true_value, col = "black", lty = 2,lwd=3) 

legend("topright", legend = "Real Effect", col = "black", lty = 2, box.lwd = 0,lwd=3)
dev.off()

m = length(names)
df = data.frame(RMSE=rep(0,m), Bias=rep(0,m), SD=rep(0,m), SE=rep(0,m), CP=rep(0,m))

for (i in 1:length(names)){
  df$Bias[i] = abs(mean(unlist(data_list[i]))-true_value)
  print(paste("absolute error: ",names[i],": ",df$Bias[i],sep = ''))
}


for (i in 1:length(names)){
  df$SD[i] = sd(unlist(data_list[i]))
  print(paste("sd:",names[i],": ",df$SD[i],sep = ''))
}

for (i in 1:length(names)){
  df$SE[i] = mean(unlist(sd_data_list[i]))
  print(paste("estimated sd:",names[i],": ",df$SE[i],sep = ''))
}

for (i in 1:length(names)){
  df$CP[i] = calculate_CP(true_value, unlist(data_list[i]), unlist(sd_data_list[i]))
  print(paste("coverage probability:",names[i],": ",df$CP[i], sep = ''))
}


for (i in 1:length(names)){
  df$RMSE[i] = mean((unlist(data_list[i])-true_value)^2)^0.5
  print(paste("RMSE: ",names[i],": ",df$RMSE[i],sep = ''))
}
df$Bias <-round(df$Bias,3)
return(df)
}



Model = "Model3"
final_df_list <- list()
randomization_methods = c("SRS","SBR","PS")
sample_sizes = c(500,1000,2000)
# randomization_methods = c("SBR")
# sample_sizes = c(500)
num_of_methods = length(randomization_methods)
num_of_samplesize = length(sample_sizes)
for (i in 1:num_of_samplesize){
  n = sample_sizes[i]
  df_list <- list()
  for (j in 1:num_of_methods) {
    df = run_simulation(Model, n, randomization_methods[j])
    df_list[[j]] <- df
  }
  df <- bind_cols(df_list)
  final_df_list[[i]] <- df
  
}

final_df = bind_rows(final_df_list)
print(final_df, row.names = FALSE)
file_name = sprintf("20251008_%s.csv", Model)
write.csv(final_df, file_name, row.names = FALSE)


