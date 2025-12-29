#! /usr/bin/env Rscript
setwd("~/semi/202311/high dimension")
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
set.seed(202311)
mc.reset.stream()


lambda = 0.75
weight = 1
omega = c(0.1,0.1,0.4,0.4)
sigma1 = 3
sigma0 = 1
pi = 1/2
Iternum = 1000
#name_methods = c("rf","nn","rpart","lasso","enet","gbm","ensemble","best")


#equal allocation
#500-50
n = 500
p = 50
name = "500_50_12"
# Model 5
Model = "Model5"
alpha0 = 1
alpha1 = 4
betavec0 = c(75,35,125,80)
betavec1 = c(100,80,60,40)
re5 = mclapply(1:Iternum,sim_DML,mc.cores = 26)
tv5 = trueval(Model,alpha0,alpha1,betavec0,betavec1)
save(re5,tv5,file = paste0("re5_",name,".RData"))

# Model 6
Model = "Model6"
alpha0 = -3
alpha1 = 0
betavec0 = c(10,24,15,20)
betavec1 = c(20,27,10)
re6 = mclapply(1:Iternum,sim_DML,mc.cores = 26)
tv6 = trueval(Model,alpha0,alpha1,betavec0,betavec1)
save(re6,tv6,file = paste0("re6_",name,".RData"))

#Model 7
Model = "Model7"
alpha0 = 5
alpha1 = 2
betavec0 = c(42,83)
betavec1 = c(30,75)
re7 = mclapply(1:Iternum,sim_DML,mc.cores = 26)
tv7 = trueval(Model,alpha0,alpha1,betavec0,betavec1)
save(re7,tv7,file = paste0("re7_",name,".RData"))

# Model 8
Model = "Model8"
alpha0 = 5
alpha1 = 5
betavec0 = c(20,30,50)
betavec1 = c(20,30,65)
re8 = mclapply(1:Iternum,sim_DML,mc.cores = 26)
tv8 = trueval(Model,alpha0,alpha1,betavec0,betavec1)
save(re8,tv8,file = paste0("re8_",name,".RData"))


#500-100
n = 500
p = 100
name = "500_100_12"
# Model 5
Model = "Model5"
alpha0 = 1
alpha1 = 4
betavec0 = c(75,35,125,80)
betavec1 = c(100,80,60,40)
re5 = mclapply(1:Iternum,sim_DML,mc.cores = 26)
tv5 = trueval(Model,alpha0,alpha1,betavec0,betavec1)
save(re5,tv5,file = paste0("re5_",name,".RData"))

# Model 6
Model = "Model6"
alpha0 = -3
alpha1 = 0
betavec0 = c(10,24,15,20)
betavec1 = c(20,27,10)
re6 = mclapply(1:Iternum,sim_DML,mc.cores = 26)
tv6 = trueval(Model,alpha0,alpha1,betavec0,betavec1)
save(re6,tv6,file = paste0("re6_",name,".RData"))

#Model 7
Model = "Model7"
alpha0 = 5
alpha1 = 2
betavec0 = c(42,83)
betavec1 = c(30,75)
re7 = mclapply(1:Iternum,sim_DML,mc.cores = 26)
tv7 = trueval(Model,alpha0,alpha1,betavec0,betavec1)
save(re7,tv7,file = paste0("re7_",name,".RData"))

# Model 8
Model = "Model8"
alpha0 = 5
alpha1 = 5
betavec0 = c(20,30,50)
betavec1 = c(20,30,65)
re8 = mclapply(1:Iternum,sim_DML,mc.cores = 26)
tv8 = trueval(Model,alpha0,alpha1,betavec0,betavec1)
save(re8,tv8,file = paste0("re8_",name,".RData"))


#500-200
n = 500
p = 200
name = "500_200_12"
# Model 5
Model = "Model5"
alpha0 = 1
alpha1 = 4
betavec0 = c(75,35,125,80)
betavec1 = c(100,80,60,40)
re5 = mclapply(1:Iternum,sim_DML,mc.cores = 26)
tv5 = trueval(Model,alpha0,alpha1,betavec0,betavec1)
save(re5,tv5,file = paste0("re5_",name,".RData"))

# Model 6
Model = "Model6"
alpha0 = -3
alpha1 = 0
betavec0 = c(10,24,15,20)
betavec1 = c(20,27,10)
re6 = mclapply(1:Iternum,sim_DML,mc.cores = 26)
tv6 = trueval(Model,alpha0,alpha1,betavec0,betavec1)
save(re6,tv6,file = paste0("re6_",name,".RData"))

#Model 7
Model = "Model7"
alpha0 = 5
alpha1 = 2
betavec0 = c(42,83)
betavec1 = c(30,75)
re7 = mclapply(1:Iternum,sim_DML,mc.cores = 26)
tv7 = trueval(Model,alpha0,alpha1,betavec0,betavec1)
save(re7,tv7,file = paste0("re7_",name,".RData"))

# Model 8
Model = "Model8"
alpha0 = 5
alpha1 = 5
betavec0 = c(20,30,50)
betavec1 = c(20,30,65)
re8 = mclapply(1:Iternum,sim_DML,mc.cores = 26)
tv8 = trueval(Model,alpha0,alpha1,betavec0,betavec1)
save(re8,tv8,file = paste0("re8_",name,".RData"))


#500-300
n = 500
p = 300
name = "500_300_12"
# Model 5
Model = "Model5"
alpha0 = 1
alpha1 = 4
betavec0 = c(75,35,125,80)
betavec1 = c(100,80,60,40)
re5 = mclapply(1:Iternum,sim_DML,mc.cores = 26)
tv5 = trueval(Model,alpha0,alpha1,betavec0,betavec1)
save(re5,tv5,file = paste0("re5_",name,".RData"))

# Model 6
Model = "Model6"
alpha0 = -3
alpha1 = 0
betavec0 = c(10,24,15,20)
betavec1 = c(20,27,10)
re6 = mclapply(1:Iternum,sim_DML,mc.cores = 26)
tv6 = trueval(Model,alpha0,alpha1,betavec0,betavec1)
save(re6,tv6,file = paste0("re6_",name,".RData"))

#Model 7
Model = "Model7"
alpha0 = 5
alpha1 = 2
betavec0 = c(42,83)
betavec1 = c(30,75)
re7 = mclapply(1:Iternum,sim_DML,mc.cores = 26)
tv7 = trueval(Model,alpha0,alpha1,betavec0,betavec1)
save(re7,tv7,file = paste0("re7_",name,".RData"))

# Model 8
Model = "Model8"
alpha0 = 5
alpha1 = 5
betavec0 = c(20,30,50)
betavec1 = c(20,30,65)
re8 = mclapply(1:Iternum,sim_DML,mc.cores = 26)
tv8 = trueval(Model,alpha0,alpha1,betavec0,betavec1)
save(re8,tv8,file = paste0("re8_",name,".RData"))


#500-500
n = 500
p = 500
name = "500_500_12"
# Model 5
Model = "Model5"
alpha0 = 1
alpha1 = 4
betavec0 = c(75,35,125,80)
betavec1 = c(100,80,60,40)
re5 = mclapply(1:Iternum,sim_DML,mc.cores = 26)
tv5 = trueval(Model,alpha0,alpha1,betavec0,betavec1)
save(re5,tv5,file = paste0("re5_",name,".RData"))

# Model 6
Model = "Model6"
alpha0 = -3
alpha1 = 0
betavec0 = c(10,24,15,20)
betavec1 = c(20,27,10)
re6 = mclapply(1:Iternum,sim_DML,mc.cores = 26)
tv6 = trueval(Model,alpha0,alpha1,betavec0,betavec1)
save(re6,tv6,file = paste0("re6_",name,".RData"))

#Model 7
Model = "Model7"
alpha0 = 5
alpha1 = 2
betavec0 = c(42,83)
betavec1 = c(30,75)
re7 = mclapply(1:Iternum,sim_DML,mc.cores = 26)
tv7 = trueval(Model,alpha0,alpha1,betavec0,betavec1)
save(re7,tv7,file = paste0("re7_",name,".RData"))

# Model 8
Model = "Model8"
alpha0 = 5
alpha1 = 5
betavec0 = c(20,30,50)
betavec1 = c(20,30,65)
re8 = mclapply(1:Iternum,sim_DML,mc.cores = 26)
tv8 = trueval(Model,alpha0,alpha1,betavec0,betavec1)
save(re8,tv8,file = paste0("re8_",name,".RData"))


#1000-50
n = 1000
p = 50
name = "1000_50_12"
# Model 5
Model = "Model5"
alpha0 = 1
alpha1 = 4
betavec0 = c(75,35,125,80)
betavec1 = c(100,80,60,40)
re5 = mclapply(1:Iternum,sim_DML,mc.cores = 26)
tv5 = trueval(Model,alpha0,alpha1,betavec0,betavec1)
save(re5,tv5,file = paste0("re5_",name,".RData"))

# Model 6
Model = "Model6"
alpha0 = -3
alpha1 = 0
betavec0 = c(10,24,15,20)
betavec1 = c(20,27,10)
re6 = mclapply(1:Iternum,sim_DML,mc.cores = 26)
tv6 = trueval(Model,alpha0,alpha1,betavec0,betavec1)
save(re6,tv6,file = paste0("re6_",name,".RData"))

#Model 7
Model = "Model7"
alpha0 = 5
alpha1 = 2
betavec0 = c(42,83)
betavec1 = c(30,75)
re7 = mclapply(1:Iternum,sim_DML,mc.cores = 26)
tv7 = trueval(Model,alpha0,alpha1,betavec0,betavec1)
save(re7,tv7,file = paste0("re7_",name,".RData"))

# Model 8
Model = "Model8"
alpha0 = 5
alpha1 = 5
betavec0 = c(20,30,50)
betavec1 = c(20,30,65)
re8 = mclapply(1:Iternum,sim_DML,mc.cores = 26)
tv8 = trueval(Model,alpha0,alpha1,betavec0,betavec1)
save(re8,tv8,file = paste0("re8_",name,".RData"))


#1000-100
n = 1000
p = 100
name = "1000_100_12"
# Model 5
Model = "Model5"
alpha0 = 1
alpha1 = 4
betavec0 = c(75,35,125,80)
betavec1 = c(100,80,60,40)
re5 = mclapply(1:Iternum,sim_DML,mc.cores = 26)
tv5 = trueval(Model,alpha0,alpha1,betavec0,betavec1)
save(re5,tv5,file = paste0("re5_",name,".RData"))

# Model 6
Model = "Model6"
alpha0 = -3
alpha1 = 0
betavec0 = c(10,24,15,20)
betavec1 = c(20,27,10)
re6 = mclapply(1:Iternum,sim_DML,mc.cores = 26)
tv6 = trueval(Model,alpha0,alpha1,betavec0,betavec1)
save(re6,tv6,file = paste0("re6_",name,".RData"))

#Model 7
Model = "Model7"
alpha0 = 5
alpha1 = 2
betavec0 = c(42,83)
betavec1 = c(30,75)
re7 = mclapply(1:Iternum,sim_DML,mc.cores = 26)
tv7 = trueval(Model,alpha0,alpha1,betavec0,betavec1)
save(re7,tv7,file = paste0("re7_",name,".RData"))

# Model 8
Model = "Model8"
alpha0 = 5
alpha1 = 5
betavec0 = c(20,30,50)
betavec1 = c(20,30,65)
re8 = mclapply(1:Iternum,sim_DML,mc.cores = 26)
tv8 = trueval(Model,alpha0,alpha1,betavec0,betavec1)
save(re8,tv8,file = paste0("re8_",name,".RData"))


#1000-200
n = 1000
p = 200
name = "1000_200_12"
# Model 5
Model = "Model5"
alpha0 = 1
alpha1 = 4
betavec0 = c(75,35,125,80)
betavec1 = c(100,80,60,40)
re5 = mclapply(1:Iternum,sim_DML,mc.cores = 26)
tv5 = trueval(Model,alpha0,alpha1,betavec0,betavec1)
save(re5,tv5,file = paste0("re5_",name,".RData"))

# Model 6
Model = "Model6"
alpha0 = -3
alpha1 = 0
betavec0 = c(10,24,15,20)
betavec1 = c(20,27,10)
re6 = mclapply(1:Iternum,sim_DML,mc.cores = 26)
tv6 = trueval(Model,alpha0,alpha1,betavec0,betavec1)
save(re6,tv6,file = paste0("re6_",name,".RData"))

#Model 7
Model = "Model7"
alpha0 = 5
alpha1 = 2
betavec0 = c(42,83)
betavec1 = c(30,75)
re7 = mclapply(1:Iternum,sim_DML,mc.cores = 26)
tv7 = trueval(Model,alpha0,alpha1,betavec0,betavec1)
save(re7,tv7,file = paste0("re7_",name,".RData"))

# Model 8
Model = "Model8"
alpha0 = 5
alpha1 = 5
betavec0 = c(20,30,50)
betavec1 = c(20,30,65)
re8 = mclapply(1:Iternum,sim_DML,mc.cores = 26)
tv8 = trueval(Model,alpha0,alpha1,betavec0,betavec1)
save(re8,tv8,file = paste0("re8_",name,".RData"))


#1000-300
n = 1000
p = 300
name = "1000_300_12"
# Model 5
Model = "Model5"
alpha0 = 1
alpha1 = 4
betavec0 = c(75,35,125,80)
betavec1 = c(100,80,60,40)
re5 = mclapply(1:Iternum,sim_DML,mc.cores = 26)
tv5 = trueval(Model,alpha0,alpha1,betavec0,betavec1)
save(re5,tv5,file = paste0("re5_",name,".RData"))

# Model 6
Model = "Model6"
alpha0 = -3
alpha1 = 0
betavec0 = c(10,24,15,20)
betavec1 = c(20,27,10)
re6 = mclapply(1:Iternum,sim_DML,mc.cores = 26)
tv6 = trueval(Model,alpha0,alpha1,betavec0,betavec1)
save(re6,tv6,file = paste0("re6_",name,".RData"))

#Model 7
Model = "Model7"
alpha0 = 5
alpha1 = 2
betavec0 = c(42,83)
betavec1 = c(30,75)
re7 = mclapply(1:Iternum,sim_DML,mc.cores = 26)
tv7 = trueval(Model,alpha0,alpha1,betavec0,betavec1)
save(re7,tv7,file = paste0("re7_",name,".RData"))

# Model 8
Model = "Model8"
alpha0 = 5
alpha1 = 5
betavec0 = c(20,30,50)
betavec1 = c(20,30,65)
re8 = mclapply(1:Iternum,sim_DML,mc.cores = 26)
tv8 = trueval(Model,alpha0,alpha1,betavec0,betavec1)
save(re8,tv8,file = paste0("re8_",name,".RData"))


#1000-500
n = 1000
p = 500
name = "1000_500_12"
# Model 5
Model = "Model5"
alpha0 = 1
alpha1 = 4
betavec0 = c(75,35,125,80)
betavec1 = c(100,80,60,40)
re5 = mclapply(1:Iternum,sim_DML,mc.cores = 26)
tv5 = trueval(Model,alpha0,alpha1,betavec0,betavec1)
save(re5,tv5,file = paste0("re5_",name,".RData"))

# Model 6
Model = "Model6"
alpha0 = -3
alpha1 = 0
betavec0 = c(10,24,15,20)
betavec1 = c(20,27,10)
re6 = mclapply(1:Iternum,sim_DML,mc.cores = 26)
tv6 = trueval(Model,alpha0,alpha1,betavec0,betavec1)
save(re6,tv6,file = paste0("re6_",name,".RData"))

#Model 7
Model = "Model7"
alpha0 = 5
alpha1 = 2
betavec0 = c(42,83)
betavec1 = c(30,75)
re7 = mclapply(1:Iternum,sim_DML,mc.cores = 26)
tv7 = trueval(Model,alpha0,alpha1,betavec0,betavec1)
save(re7,tv7,file = paste0("re7_",name,".RData"))

# Model 8
Model = "Model8"
alpha0 = 5
alpha1 = 5
betavec0 = c(20,30,50)
betavec1 = c(20,30,65)
re8 = mclapply(1:Iternum,sim_DML,mc.cores = 26)
tv8 = trueval(Model,alpha0,alpha1,betavec0,betavec1)
save(re8,tv8,file = paste0("re8_",name,".RData"))


#unequal allocation
pi = 2/3

#500-50
n = 500
p = 50
name = "500_50_23"
# Model 5
Model = "Model5"
alpha0 = 1
alpha1 = 4
betavec0 = c(75,35,125,80)
betavec1 = c(100,80,60,40)
re5 = mclapply(1:Iternum,sim_DML,mc.cores = 26)
tv5 = trueval(Model,alpha0,alpha1,betavec0,betavec1)
save(re5,tv5,file = paste0("re5_",name,".RData"))

# Model 6
Model = "Model6"
alpha0 = -3
alpha1 = 0
betavec0 = c(10,24,15,20)
betavec1 = c(20,27,10)
re6 = mclapply(1:Iternum,sim_DML,mc.cores = 26)
tv6 = trueval(Model,alpha0,alpha1,betavec0,betavec1)
save(re6,tv6,file = paste0("re6_",name,".RData"))

#Model 7
Model = "Model7"
alpha0 = 5
alpha1 = 2
betavec0 = c(42,83)
betavec1 = c(30,75)
re7 = mclapply(1:Iternum,sim_DML,mc.cores = 26)
tv7 = trueval(Model,alpha0,alpha1,betavec0,betavec1)
save(re7,tv7,file = paste0("re7_",name,".RData"))

# Model 8
Model = "Model8"
alpha0 = 5
alpha1 = 5
betavec0 = c(20,30,50)
betavec1 = c(20,30,65)
re8 = mclapply(1:Iternum,sim_DML,mc.cores = 26)
tv8 = trueval(Model,alpha0,alpha1,betavec0,betavec1)
save(re8,tv8,file = paste0("re8_",name,".RData"))


#500-100
n = 500
p = 100
name = "500_100_23"

# Model 5
Model = "Model5"
alpha0 = 1
alpha1 = 4
betavec0 = c(75,35,125,80)
betavec1 = c(100,80,60,40)
re5 = mclapply(1:Iternum,sim_DML,mc.cores = 26)
tv5 = trueval(Model,alpha0,alpha1,betavec0,betavec1)
save(re5,tv5,file = paste0("re5_",name,".RData"))

# Model 6
Model = "Model6"
alpha0 = -3
alpha1 = 0
betavec0 = c(10,24,15,20)
betavec1 = c(20,27,10)
re6 = mclapply(1:Iternum,sim_DML,mc.cores = 26)
tv6 = trueval(Model,alpha0,alpha1,betavec0,betavec1)
save(re6,tv6,file = paste0("re6_",name,".RData"))

#Model 7
Model = "Model7"
alpha0 = 5
alpha1 = 2
betavec0 = c(42,83)
betavec1 = c(30,75)
re7 = mclapply(1:Iternum,sim_DML,mc.cores = 26)
tv7 = trueval(Model,alpha0,alpha1,betavec0,betavec1)
save(re7,tv7,file = paste0("re7_",name,".RData"))

# Model 8
Model = "Model8"
alpha0 = 5
alpha1 = 5
betavec0 = c(20,30,50)
betavec1 = c(20,30,65)
re8 = mclapply(1:Iternum,sim_DML,mc.cores = 26)
tv8 = trueval(Model,alpha0,alpha1,betavec0,betavec1)
save(re8,tv8,file = paste0("re8_",name,".RData"))


#500-200
n = 500
p = 200
name = "500_200_23"

# Model 5
Model = "Model5"
alpha0 = 1
alpha1 = 4
betavec0 = c(75,35,125,80)
betavec1 = c(100,80,60,40)
re5 = mclapply(1:Iternum,sim_DML,mc.cores = 26)
tv5 = trueval(Model,alpha0,alpha1,betavec0,betavec1)
save(re5,tv5,file = paste0("re5_",name,".RData"))

# Model 6
Model = "Model6"
alpha0 = -3
alpha1 = 0
betavec0 = c(10,24,15,20)
betavec1 = c(20,27,10)
re6 = mclapply(1:Iternum,sim_DML,mc.cores = 26)
tv6 = trueval(Model,alpha0,alpha1,betavec0,betavec1)
save(re6,tv6,file = paste0("re6_",name,".RData"))

#Model 7
Model = "Model7"
alpha0 = 5
alpha1 = 2
betavec0 = c(42,83)
betavec1 = c(30,75)
re7 = mclapply(1:Iternum,sim_DML,mc.cores = 26)
tv7 = trueval(Model,alpha0,alpha1,betavec0,betavec1)
save(re7,tv7,file = paste0("re7_",name,".RData"))

# Model 8
Model = "Model8"
alpha0 = 5
alpha1 = 5
betavec0 = c(20,30,50)
betavec1 = c(20,30,65)
re8 = mclapply(1:Iternum,sim_DML,mc.cores = 26)
tv8 = trueval(Model,alpha0,alpha1,betavec0,betavec1)
save(re8,tv8,file = paste0("re8_",name,".RData"))

#500-300
n = 500
p = 300
name = "500_300_23"

# Model 5
Model = "Model5"
alpha0 = 1
alpha1 = 4
betavec0 = c(75,35,125,80)
betavec1 = c(100,80,60,40)
re5 = mclapply(1:Iternum,sim_DML,mc.cores = 26)
tv5 = trueval(Model,alpha0,alpha1,betavec0,betavec1)
save(re5,tv5,file = paste0("re5_",name,".RData"))

# Model 6
Model = "Model6"
alpha0 = -3
alpha1 = 0
betavec0 = c(10,24,15,20)
betavec1 = c(20,27,10)
re6 = mclapply(1:Iternum,sim_DML,mc.cores = 26)
tv6 = trueval(Model,alpha0,alpha1,betavec0,betavec1)
save(re6,tv6,file = paste0("re6_",name,".RData"))

#Model 7
Model = "Model7"
alpha0 = 5
alpha1 = 2
betavec0 = c(42,83)
betavec1 = c(30,75)
re7 = mclapply(1:Iternum,sim_DML,mc.cores = 26)
tv7 = trueval(Model,alpha0,alpha1,betavec0,betavec1)
save(re7,tv7,file = paste0("re7_",name,".RData"))

# Model 8
Model = "Model8"
alpha0 = 5
alpha1 = 5
betavec0 = c(20,30,50)
betavec1 = c(20,30,65)
re8 = mclapply(1:Iternum,sim_DML,mc.cores = 26)
tv8 = trueval(Model,alpha0,alpha1,betavec0,betavec1)
save(re8,tv8,file = paste0("re8_",name,".RData"))

#500-500
n = 500
p = 500
name = "500_500_23"

# Model 5
Model = "Model5"
alpha0 = 1
alpha1 = 4
betavec0 = c(75,35,125,80)
betavec1 = c(100,80,60,40)
re5 = mclapply(1:Iternum,sim_DML,mc.cores = 26)
tv5 = trueval(Model,alpha0,alpha1,betavec0,betavec1)
save(re5,tv5,file = paste0("re5_",name,".RData"))

# Model 6
Model = "Model6"
alpha0 = -3
alpha1 = 0
betavec0 = c(10,24,15,20)
betavec1 = c(20,27,10)
re6 = mclapply(1:Iternum,sim_DML,mc.cores = 26)
tv6 = trueval(Model,alpha0,alpha1,betavec0,betavec1)
save(re6,tv6,file = paste0("re6_",name,".RData"))

#Model 7
Model = "Model7"
alpha0 = 5
alpha1 = 2
betavec0 = c(42,83)
betavec1 = c(30,75)
re7 = mclapply(1:Iternum,sim_DML,mc.cores = 26)
tv7 = trueval(Model,alpha0,alpha1,betavec0,betavec1)
save(re7,tv7,file = paste0("re7_",name,".RData"))

# Model 8
Model = "Model8"
alpha0 = 5
alpha1 = 5
betavec0 = c(20,30,50)
betavec1 = c(20,30,65)
re8 = mclapply(1:Iternum,sim_DML,mc.cores = 26)
tv8 = trueval(Model,alpha0,alpha1,betavec0,betavec1)
save(re8,tv8,file = paste0("re8_",name,".RData"))

#1000-50
n = 1000
p = 50
name = "1000_50_23"

# Model 5
Model = "Model5"
alpha0 = 1
alpha1 = 4
betavec0 = c(75,35,125,80)
betavec1 = c(100,80,60,40)
re5 = mclapply(1:Iternum,sim_DML,mc.cores = 26)
tv5 = trueval(Model,alpha0,alpha1,betavec0,betavec1)
save(re5,tv5,file = paste0("re5_",name,".RData"))

# Model 6
Model = "Model6"
alpha0 = -3
alpha1 = 0
betavec0 = c(10,24,15,20)
betavec1 = c(20,27,10)
re6 = mclapply(1:Iternum,sim_DML,mc.cores = 26)
tv6 = trueval(Model,alpha0,alpha1,betavec0,betavec1)
save(re6,tv6,file = paste0("re6_",name,".RData"))

#Model 7
Model = "Model7"
alpha0 = 5
alpha1 = 2
betavec0 = c(42,83)
betavec1 = c(30,75)
re7 = mclapply(1:Iternum,sim_DML,mc.cores = 26)
tv7 = trueval(Model,alpha0,alpha1,betavec0,betavec1)
save(re7,tv7,file = paste0("re7_",name,".RData"))

# Model 8
Model = "Model8"
alpha0 = 5
alpha1 = 5
betavec0 = c(20,30,50)
betavec1 = c(20,30,65)
re8 = mclapply(1:Iternum,sim_DML,mc.cores = 26)
tv8 = trueval(Model,alpha0,alpha1,betavec0,betavec1)
save(re8,tv8,file = paste0("re8_",name,".RData"))

#1000-100
n = 1000
p = 100
name = "1000_100_23"

# Model 5
Model = "Model5"
alpha0 = 1
alpha1 = 4
betavec0 = c(75,35,125,80)
betavec1 = c(100,80,60,40)
re5 = mclapply(1:Iternum,sim_DML,mc.cores = 26)
tv5 = trueval(Model,alpha0,alpha1,betavec0,betavec1)
save(re5,tv5,file = paste0("re5_",name,".RData"))

# Model 6
Model = "Model6"
alpha0 = -3
alpha1 = 0
betavec0 = c(10,24,15,20)
betavec1 = c(20,27,10)
re6 = mclapply(1:Iternum,sim_DML,mc.cores = 26)
tv6 = trueval(Model,alpha0,alpha1,betavec0,betavec1)
save(re6,tv6,file = paste0("re6_",name,".RData"))

#Model 7
Model = "Model7"
alpha0 = 5
alpha1 = 2
betavec0 = c(42,83)
betavec1 = c(30,75)
re7 = mclapply(1:Iternum,sim_DML,mc.cores = 26)
tv7 = trueval(Model,alpha0,alpha1,betavec0,betavec1)
save(re7,tv7,file = paste0("re7_",name,".RData"))

# Model 8
Model = "Model8"
alpha0 = 5
alpha1 = 5
betavec0 = c(20,30,50)
betavec1 = c(20,30,65)
re8 = mclapply(1:Iternum,sim_DML,mc.cores = 26)
tv8 = trueval(Model,alpha0,alpha1,betavec0,betavec1)
save(re8,tv8,file = paste0("re8_",name,".RData"))

#1000-200
n = 1000
p = 200
name = "1000_200_23"

# Model 5
Model = "Model5"
alpha0 = 1
alpha1 = 4
betavec0 = c(75,35,125,80)
betavec1 = c(100,80,60,40)
re5 = mclapply(1:Iternum,sim_DML,mc.cores = 26)
tv5 = trueval(Model,alpha0,alpha1,betavec0,betavec1)
save(re5,tv5,file = paste0("re5_",name,".RData"))

# Model 6
Model = "Model6"
alpha0 = -3
alpha1 = 0
betavec0 = c(10,24,15,20)
betavec1 = c(20,27,10)
re6 = mclapply(1:Iternum,sim_DML,mc.cores = 26)
tv6 = trueval(Model,alpha0,alpha1,betavec0,betavec1)
save(re6,tv6,file = paste0("re6_",name,".RData"))

#Model 7
Model = "Model7"
alpha0 = 5
alpha1 = 2
betavec0 = c(42,83)
betavec1 = c(30,75)
re7 = mclapply(1:Iternum,sim_DML,mc.cores = 26)
tv7 = trueval(Model,alpha0,alpha1,betavec0,betavec1)
save(re7,tv7,file = paste0("re7_",name,".RData"))

# Model 8
Model = "Model8"
alpha0 = 5
alpha1 = 5
betavec0 = c(20,30,50)
betavec1 = c(20,30,65)
re8 = mclapply(1:Iternum,sim_DML,mc.cores = 26)
tv8 = trueval(Model,alpha0,alpha1,betavec0,betavec1)
save(re8,tv8,file = paste0("re8_",name,".RData"))

#1000-300
n = 1000
p = 300
name = "1000_300_23"

# Model 5
Model = "Model5"
alpha0 = 1
alpha1 = 4
betavec0 = c(75,35,125,80)
betavec1 = c(100,80,60,40)
re5 = mclapply(1:Iternum,sim_DML,mc.cores = 26)
tv5 = trueval(Model,alpha0,alpha1,betavec0,betavec1)
save(re5,tv5,file = paste0("re5_",name,".RData"))

# Model 6
Model = "Model6"
alpha0 = -3
alpha1 = 0
betavec0 = c(10,24,15,20)
betavec1 = c(20,27,10)
re6 = mclapply(1:Iternum,sim_DML,mc.cores = 26)
tv6 = trueval(Model,alpha0,alpha1,betavec0,betavec1)
save(re6,tv6,file = paste0("re6_",name,".RData"))

#Model 7
Model = "Model7"
alpha0 = 5
alpha1 = 2
betavec0 = c(42,83)
betavec1 = c(30,75)
re7 = mclapply(1:Iternum,sim_DML,mc.cores = 26)
tv7 = trueval(Model,alpha0,alpha1,betavec0,betavec1)
save(re7,tv7,file = paste0("re7_",name,".RData"))

# Model 8
Model = "Model8"
alpha0 = 5
alpha1 = 5
betavec0 = c(20,30,50)
betavec1 = c(20,30,65)
re8 = mclapply(1:Iternum,sim_DML,mc.cores = 26)
tv8 = trueval(Model,alpha0,alpha1,betavec0,betavec1)
save(re8,tv8,file = paste0("re8_",name,".RData"))

#1000-500
n = 1000
p = 500
name = "1000_500_23"
# Model 5
Model = "Model5"
alpha0 = 1
alpha1 = 4
betavec0 = c(75,35,125,80)
betavec1 = c(100,80,60,40)
re5 = mclapply(1:Iternum,sim_DML,mc.cores = 26)
tv5 = trueval(Model,alpha0,alpha1,betavec0,betavec1)
save(re5,tv5,file = paste0("re5_",name,".RData"))

# Model 6
Model = "Model6"
alpha0 = -3
alpha1 = 0
betavec0 = c(10,24,15,20)
betavec1 = c(20,27,10)
re6 = mclapply(1:Iternum,sim_DML,mc.cores = 26)
tv6 = trueval(Model,alpha0,alpha1,betavec0,betavec1)
save(re6,tv6,file = paste0("re6_",name,".RData"))

#Model 7
Model = "Model7"
alpha0 = 5
alpha1 = 2
betavec0 = c(42,83)
betavec1 = c(30,75)
re7 = mclapply(1:Iternum,sim_DML,mc.cores = 26)
tv7 = trueval(Model,alpha0,alpha1,betavec0,betavec1)
save(re7,tv7,file = paste0("re7_",name,".RData"))

# Model 8
Model = "Model8"
alpha0 = 5
alpha1 = 5
betavec0 = c(20,30,50)
betavec1 = c(20,30,65)
re8 = mclapply(1:Iternum,sim_DML,mc.cores = 26)
tv8 = trueval(Model,alpha0,alpha1,betavec0,betavec1)
save(re8,tv8,file = paste0("re8_",name,".RData"))
