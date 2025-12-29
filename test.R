#! /usr/bin/env Rscript
#setwd("~/semi/202209/0907")
library(microbenchmark)
Rcpp::sourceCpp("vestimator.cpp")
source("estimators.R")
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
set.seed(202209)
#mc.reset.stream()

n = 400
lambda = 0.75
weight = 1
omega = c(0.1,0.1,0.4,0.4)
sigma1 = 3
sigma0 = 1
pi = 1/2
Iternum = 1
p = 200
#name_methods = c("rf","nn","rpart","lasso","enet","gbm","ensemble","best")

# Model 5
Model = "Model5"
alpha0 = 1
alpha1 = 4
betavec0 = c(75,35,125,80)
betavec1 = c(100,80,60,40)
result = sim_DML(1)
