source("cal_solve_weights.R")
library(MASS)


cal_est<- function(strt_num,str_indicator,A,pai_each_data,hX_1_str,hX_0_str,Y,agg=FALSE,standard_weights=TRUE,EL=FALSE,uniform_weights=FALSE,population_pi=NA){
  
  
  if (uniform_weights){
    cal_weights = rep(1,length(Y))
  }
  else {
    cal_weights = cal_solve_weights_separate(strt_num,str_indicator,A,pai_each_data,hX_1_str,hX_0_str,EL=EL)
    #print(cal_weights)
  }
  
  
  Y_str_mean = Y
  cal_weights_mean = cal_weights
  
  for(i in 1:strt_num){
    ind1 = which(str_indicator == i & A == 1)
    ind0 = which(str_indicator == i & A == 0)
    ind_all = which(str_indicator == i)
    Y_str_mean[ind0] = mean(Y[ind0])
    Y_str_mean[ind1] = mean(Y[ind1])
    cal_weights_mean[ind_all] = mean(cal_weights[ind_all])
    
  }
  
  
  interaction_est = mean(((A / pai_each_data)-((1-A)/(1-pai_each_data))) * Y)
  adjust_term = mean(cal_weights * ((A / pai_each_data)-((1-A)/(1-pai_each_data))) * (Y-Y_str_mean))
  tau_hat_str = interaction_est + adjust_term
  # tau_hat_str = mean(cal_weights * ((A / pai_each_data)-((1-A)/(1-pai_each_data))) * (Y))
  
  
  
  # standard error
  dof = ncol(hX_1_str) + ncol(hX_0_str) + 1 - all(hX_1_str == 0) - all(hX_0_str == 0)
  n = length(Y)
  Ybar1 = mean(Y[which(A == 1)])
  Ybar0 = mean(Y[which(A == 0)])
  hx = cbind(hX_1_str,hX_0_str)
  v = 0
  for(i in 1:strt_num){
    ind1 = which(str_indicator == i & A == 1)
    ind0 = which(str_indicator == i & A == 0)
    ind_all = which(str_indicator == i)
    
    nk = length(ind_all)
    n1k = length(ind1)
    n0k = length(ind0)
    pink = n1k / nk
    if ((length(ind1) <= 3 | length(ind0) <= 3) & !is.na(population_pi)){
      
      pink = population_pi
      }
    
    Ybar0k = mean(Y[ind0])
    Ybar1k = mean(Y[ind1])
    
    pnk =  nk / n
    v1 = pnk / (1-pink) * var(Y[ind0]) + pnk / pink * var(Y[ind1])
    v2 = pnk * (Ybar1k-Ybar0k-Ybar1+Ybar0)^2
    
    if (uniform_weights | (all(cal_weights == 1))) {
      v = v + (v1 + v2)
    }
    else {
      adjusted_factor = nk / (nk - dof)
      hxk_demean = scale(hx[ind_all,], center = TRUE, scale = FALSE)
      tranformed_y = (1-pink) / pink * A[ind_all] * (Y[ind_all] - Ybar1k) +  pink / (1-pink) * (1-A[ind_all]) * (Y[ind_all] - Ybar0k)
      XTY = t(hxk_demean) %*% tranformed_y / n
      XTX = t(hxk_demean) %*% diag((A[ind_all]-pink)^2) %*% hxk_demean / n
      v3 = t(XTY) %*% ginv(XTX) %*% XTY
      v = v + adjusted_factor * (v1 + v2 -v3)
    }
    
  }
  
  
  return(c(tau_hat_str,sqrt(v/n)))
}



