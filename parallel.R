library(parallel)

sim_fun<-function(Model,method){
  return(match.fun(paste(Model,"_",method,sep = '')))
}

sim_low.in <- function(data,Model){
  ATE <- matrix(0,nrow=2,ncol=5)
  sd <- matrix(0,nrow=2,ncol=5)
  result_lin <- est_lin(data)
  result_kernel <- est_kernel(data)
  result_ns <- est_ns(data,Model)
  result_nnet <- est_nnet(data)
  result_rf <- est_rf(data)
  ATE[,1] <- result_lin[,1]
  ATE[,2] <- result_kernel[,1]
  ATE[,3] <- result_ns[,1]
  ATE[,4] <- result_nnet[,1]
  ATE[,5] <- result_rf[,1]
  sd[,1] <- result_lin[,2]
  sd[,2] <- result_kernel[,2]
  sd[,3] <- result_ns[,2]
  sd[,4] <- result_nnet[,2]
  sd[,5] <- result_rf[,2]
  rownames(ATE) <- c("common","specific")
  colnames(ATE) <- c("linear","kernel","ns","nnet","rf")
  rownames(sd) <- c("common","specific")
  colnames(sd) <- c("linear","kernel","ns","nnet","rf")
  return(list(ATE = ATE,sd = sd))
}

sim_low <- function(it){
  FUN = sim_fun(Model,"SRS")
  data = FUN(n,alpha0,alpha1,betavec0,betavec1,pi,sigma0,sigma1)
  result_SRS <- sim_low.in(data,Model)
  
  FUN = sim_fun(Model,"SBR")
  data = FUN(n,alpha0,alpha1,betavec0,betavec1,pi,sigma0,sigma1)
  result_SBR <- sim_low.in(data,Model)
  
  FUN = sim_fun(Model,"PS")
  data = FUN(n,alpha0,alpha1,betavec0,betavec1,pi,sigma0,sigma1,weight,lambda)
  result_PS <- sim_low.in(data,Model)
  
  return(list(result_SRS = result_SRS,result_SBR = result_SBR,result_PS = result_PS))
}

sim_high.in <- function(data,Model,hidden){
  ATE <- matrix(0,nrow=2,ncol=5)
  sd <- matrix(0,nrow=2,ncol=5)
  result_lasso <- est_lasso(data)
  result_nnet <- est_neural(data,2)
  #result_nnet_deep <- est_neural(data,rep(200,10))
  result_gbm <- est_gbm(data)
  result_rf <- est_rf(data)
  result_rpart <- est_cart(data)
  ATE[,1] <- result_lasso[,1]
  ATE[,2] <- result_rf[,1]
  ATE[,3] <- result_gbm[,1]
  ATE[,4] <- result_rpart[,1]
  ATE[,5] <- result_nnet[,1]
  #ATE[,6] <- result_nnet_deep[,1]
  sd[,1] <- result_lasso[,2]
  sd[,2] <- result_rf[,2]
  sd[,3] <- result_gbm[,2]
  sd[,4] <- result_rpart[,2]
  sd[,5] <- result_nnet[,2]
  #sd[,6] <- result_nnet_deep[,2]
  rownames(ATE) <- c("common","specific")
  colnames(ATE) <- c("lasso","rf","gbm","rpart","neural net")
  rownames(sd) <- c("common","specific")
  colnames(sd) <- c("lasso","rf","gbm","rpart","neural net")
  return(list(ATE = ATE,sd = sd))
}

sim_high <- function(it){
  FUN = sim_fun(Model,"SRS")
  data = FUN(n,alpha0,alpha1,betavec0,betavec1,pi,sigma0,sigma1,p)
  result_SRS <- sim_high.in(data,Model)
  
  FUN = sim_fun(Model,"SBR")
  data = FUN(n,alpha0,alpha1,betavec0,betavec1,pi,sigma0,sigma1,p)
  result_SBR <- sim_high.in(data,Model)
  
  FUN = sim_fun(Model,"PS")
  data = FUN(n,alpha0,alpha1,betavec0,betavec1,pi,sigma0,sigma1,p,weight,lambda)
  result_PS <- sim_high.in(data,Model)
  
  return(list(result_SRS = result_SRS,result_SBR = result_SBR,result_PS = result_PS))
}


sim_DML.in <- function(data,Model){
  ATE <- matrix(0,nrow=2,ncol=5)
  sd <- matrix(0,nrow=2,ncol=5)
  result_lasso <- DML_lasso(data)
  result_nnet <- DML_neural(data,2)
  #result_nnet_deep <- DML_neural(data,rep(200,10))
  result_gbm <- DML_gbm(data)
  result_rf <- DML_rf(data)
  result_rpart <- DML_cart(data)
  ATE[,1] <- result_lasso[,1]
  ATE[,2] <- result_rf[,1]
  ATE[,3] <- result_gbm[,1]
  ATE[,4] <- result_rpart[,1]
  ATE[,5] <- result_nnet[,1]
  #ATE[,6] <- result_nnet_deep[,1]
  sd[,1] <- result_lasso[,2]
  sd[,2] <- result_rf[,2]
  sd[,3] <- result_gbm[,2]
  sd[,4] <- result_rpart[,2]
  sd[,5] <- result_nnet[,2]
  #sd[,6] <- result_nnet_deep[,2]
  rownames(ATE) <- c("common","specific")
  colnames(ATE) <- c("lasso","rf","gbm","rpart","neural net")
  rownames(sd) <- c("common","specific")
  colnames(sd) <- c("lasso","rf","gbm","rpart","neural net")
  return(list(ATE = ATE,sd = sd))
}

sim_DML <- function(it){
  FUN = sim_fun(Model,"SRS")
  data = FUN(n,alpha0,alpha1,betavec0,betavec1,pi,sigma0,sigma1,p)
  print(data)
  result_SRS <- sim_high.in(data,Model)
  result_SRS_DML <- sim_DML.in(data,Model)
  
  FUN = sim_fun(Model,"SBR")
  data = FUN(n,alpha0,alpha1,betavec0,betavec1,pi,sigma0,sigma1,p)
  result_SBR <- sim_high.in(data,Model)
  result_SBR_DML <- sim_DML.in(data,Model)
  
  FUN = sim_fun(Model,"PS")
  data = FUN(n,alpha0,alpha1,betavec0,betavec1,pi,sigma0,sigma1,p,weight,lambda)
  result_PS <- sim_high.in(data,Model)
  result_PS_DML <- sim_DML.in(data,Model)
  
  return(list(result_SRS = result_SRS,result_SBR = result_SBR,result_PS = result_PS,
              result_SRS_DML = result_SRS_DML,result_SBR_DML = result_SBR_DML,result_PS_DML = result_PS_DML))
}

sim_ensemble<-function(it){
  
  FUN = sim_fun(Model,"SRS")
  data = FUN(n,alpha0,alpha1,betavec0,betavec1,pi,sigma0,sigma1,p)
  result_SRS = DML_ensemble(data,name_methods)
  
  FUN = sim_fun(Model,"SBR")
  data = FUN(n,alpha0,alpha1,betavec0,betavec1,pi,sigma0,sigma1,p)
  result_SBR = DML_ensemble(data,name_methods)

  
  return(list(result_SRS = result_SRS,result_SBR = result_SBR))
}

cp_ensemble<-function(result,tv,name_methods){
  Iternum = length(result)
  lnm = length(name_methods)
  SRS_proposal_cp = numeric(lnm)
  SRS_strat_cp = numeric(lnm)
  SBR_proposal_cp = numeric(lnm)
  SBR_strat_cp = numeric(lnm)
  for(i in 1:Iternum){
    SRS_proposal_cp = SRS_proposal_cp + (result[[i]]$result_SRS$res_tau[1,]+1.96*result[[i]]$result_SRS$res_var[1,]>=tv)*
      (result[[i]]$result_SRS$res_tau[1,]-1.96*result[[i]]$result_SRS$res_var[1,]<=tv)
    SRS_strat_cp = SRS_strat_cp + (result[[i]]$result_SRS$res_tau[2,]+1.96*result[[i]]$result_SRS$res_var[2,]>=tv)*
      (result[[i]]$result_SRS$res_tau[2,]-1.96*result[[i]]$result_SRS$res_var[2,]<=tv)
    SBR_proposal_cp = SBR_proposal_cp + (result[[i]]$result_SBR$res_tau[1,]+1.96*result[[i]]$result_SBR$res_var[1,]>=tv)*
      (result[[i]]$result_SBR$res_tau[1,]-1.96*result[[i]]$result_SBR$res_var[1,]<=tv)
    SBR_strat_cp = SBR_strat_cp + (result[[i]]$result_SBR$res_tau[2,]+1.96*result[[i]]$result_SBR$res_var[2,]>=tv)*
      (result[[i]]$result_SBR$res_tau[2,]-1.96*result[[i]]$result_SBR$res_var[2,]<=tv)
  }
  result = rbind(SRS_proposal_cp,SRS_strat_cp,SBR_proposal_cp,SBR_strat_cp)/Iternum
  rownames(result) = c("SRS_pro","SRS_strat","SBR_pro","SBR_strat")
  colnames(result) = name_methods
  return(result)
}
