

library(pls)
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
source("Model5.R")
source("Model6.R")
source("Model7.R")
source("Model8.R")
source("tables.R")
source("DML.R")
source("parallel.R")
RNGkind("L'Ecuyer-CMRG")
set.seed(202311)
Model = "Model6"
alpha0 = 1
alpha1 = 4
betavec0 = c(75,35,125,80)
betavec1 = c(100,80,60,40)
FUN = sim_fun(Model,"SRS")
lambda = 0.75
weight = 1
omega = c(0.1,0.1,0.4,0.4)
sigma1 = 3
sigma0 = 1
pi = 0.5

p=30
S=300

n = 100000
data = FUN(n,alpha0,alpha1,betavec0,betavec1,pi,sigma0,sigma1,p)
# result = cal_rf(data)
effect = mean(2 * (1-data$A) * (data$S==data$S[1]) * data$Y)
pi1 = length(which(data$S==data$S[1] & data$A == 1))  / length(which(data$S==data$S[1]))
Y =  (1-data$A) * (data$S==data$S[1]) * data$Y / pi1- effect
strt = c(unique(data$S))
strt_num = length(strt)
X = matrix(0,nrow = n, ncol=2 * strt_num)
for (i in 1:n){
  str_for_this_sample = data$S[i]
  treat_for_this_sample = data$A[i]
  indica = which(strt == str_for_this_sample)
  X[i,indica] = ( (1-treat_for_this_sample) / (1 - pi1) - 1) * data$X[i,1]
  X[i,indica+strt_num] = (treat_for_this_sample / pi1 - 1) * data$X[i,2]
}
model_EL = lm(Y~X-1)
coef_EL = model_EL$coefficients
summary(model_EL)
# 
# Y_cal = 2 * (1-data$A) * (data$S==data$S[1]) * data$Y
# X_cal = matrix(0,nrow = n, ncol=2 * strt_num)
# for (i in 1:n){
#   str_for_this_sample = data$S[i]
#   treat_for_this_sample = data$A[i]
#   
#   ind = which(strt == str_for_this_sample)
#   X_cal[i,ind] = 2 * (1-treat_for_this_sample) * data$X[i,1]
#   X_cal[i,ind+strt_num] = 2 * treat_for_this_sample * data$X[i,2]
# }
# model_cal = lm(Y_cal~X_cal-1)
# coef_cal = model_cal$coefficients
# summary(model_cal)





X = rnorm(10000,0,1)

Y = X + rnorm(10000,0,0.2)+3
intercept = rep(1,10000)
Yii = Y-mean(Y)
Xii = X-mean(X)
model = lm(Yii~Xii-1)
summary(model)
model = lm(Y~Xii-1)
summary(model)
# 
# Y_demean = Y-mean(Y)
# X_demean =  X - colMeans(X)
# model_2 = lm(Y_demean~X_demean-1)
# summary(model_2)