library(Rcpp)
library(MASS)
sourceCpp("generate_objective_func.cpp")


cal_est<- function(strt_num,str_indicator,A,pai_each_data,hX_1_str,hX_0_str,Y,EL=FALSE,uniform_weights=FALSE,population_pi=NA){
  
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



d.rho <- function(v) {
  # return(exp(-v))
  # return(-v)
  return(1 / (v + 1))
}

rho <- function(v) {
  # return(-exp(-v))
  # return(-v^2/2)
  return(log(v + 1))
}


cal_solve_weights <- function(str_num, str_indicator, A, pai, hX_1, hX_0, EL = FALSE) {
  if (all(hX_1 == 0) & all(hX_0 == 0)) {
    return(rep(1, length(A)))
  }

  # hX_1 <- scale(hX_1)
  # hX_0 <- scale(hX_0)
  # print(hX_1)
  # print(hX_0)
  square <- TRUE
  # print("hX_1:")
  # print(hX_1)
  # print("HX1_end")

  result <- generate_objective_func(
    str_num,
    str_indicator,
    A,
    pai,
    as.matrix(hX_1),
    as.matrix(hX_0)
  )


  if (EL) {
    start_values <- rnorm(length(result$B), 0, 0.0001)

    # start_values = rep(0,length(result$B))

    opti_result <- optim(
      par = start_values,
      fn = function(x) neg_objective_func(x, result$mat_A, result$B),
      control = list(
        reltol = 1e-16,
        maxit = 10000
      )
    )

    lam_optimal <- opti_result$par
    # print(lam_optimal)
    estimated_weights <- t(d.rho(t(lam_optimal) %*% result$mat_A))
    # print(max(estimated_weights))
    # print(min(estimated_weights))
    # print((result$mat_A %*% estimated_weights) / length(A)  - result$B)
    return(estimated_weights)
  }

  if (square) {
    # print(result$mat_A)
    XTX <- result$mat_A %*% t(result$mat_A) / length(A)
    # print(result$mat_A)
    # print(result$B)
    # print(XTX)
    # print(rowMeans(result$mat_A) - result$B)
    lam <- ginv(XTX) %*% (rowMeans(result$mat_A) - result$B)
    # print(lam)
    estimated_weights <- -t(t(lam) %*% result$mat_A) + 1
    # print(t(estimated_weights))
    # print((result$mat_A %*% estimated_weights) / length(A)  - result$B)
    # print(estimated_weights)
    return(estimated_weights)
  }
}


neg_objective_func <- function(lam, mat_A, B) {
  return(-mean(rho(t(lam) %*% mat_A)) + lam %*% B)
}


cal_solve_weights_separate <- function(str_num, str_indicator, A, pai, hX_1, hX_0, EL = FALSE) {
  n <- length(A)
  weights <- rep(1, n)

  for (i in 1:str_num) {
    indistr <- which(str_indicator == i)
    hX_1 <- as.matrix(hX_1)
    hX_0 <- as.matrix(hX_0)
    weightsistr <- cal_solve_weights(1, rep(1, length(indistr)),
      A[indistr],
      pai[indistr],
      hX_1[indistr, ],
      hX_0[indistr, ],
      EL = EL
    )
    weights[indistr] <- weightsistr
  }
  return(weights)
}
