library(Rcpp)
library(MASS)
sourceCpp("generate_objective_func.cpp")

d.rho<-function(v){
  #return(exp(-v))
  #return(-v)
  return(1/(v+1))
  
}

rho<-function(v){
  #return(-exp(-v))
  # return(-v^2/2)
  return(log(v+1))
}







cal_solve_weights <- function(str_num, str_indicator, A, pai, hX_1, hX_0, standard_weights=TRUE, EL=FALSE){

  if (all(hX_1 == 0) & all(hX_0 == 0)) {
    return(rep(1,length(A)))
  }
  
  # hX_1 <- scale(hX_1)
  # hX_0 <- scale(hX_0)
  # print(hX_1)
  # print(hX_0)
  square = TRUE
  #print("hX_1:")
  #print(hX_1)
  #print("HX1_end")
  
  result <- generate_objective_func(str_num, 
                                      str_indicator, 
                                      A, 
                                      pai, 
                                      as.matrix(hX_1), 
                                      as.matrix(hX_0))
  
  
  
  
  
  if (EL) {
    
    start_values = rnorm(length(result$B),0,0.0001)
    
    # start_values = rep(0,length(result$B))
    
    opti_result = optim(par = start_values, 
                        fn = function(x) neg_objective_func(x, result$mat_A,result$B),
                        control = list(
                          reltol = 1e-16,              
                          maxit = 10000               
                        )
    )
    
    lam_optimal = opti_result$par
    #print(lam_optimal)
    estimated_weights = t(d.rho(t(lam_optimal) %*% result$mat_A))
    #print(max(estimated_weights))
    #print(min(estimated_weights))
    #print((result$mat_A %*% estimated_weights) / length(A)  - result$B)
    return(estimated_weights)
  
  }
  
  if (square) {
    #print(result$mat_A)
    XTX = result$mat_A %*% t(result$mat_A) / length(A)
    #print(result$mat_A)
    #print(result$B)
    # print(XTX)
    # print(rowMeans(result$mat_A) - result$B)
    lam = ginv(XTX)%*% (rowMeans(result$mat_A) - result$B)
    #print(lam)
    estimated_weights = -t(t(lam) %*% result$mat_A) + 1
    #print(t(estimated_weights))
    #print((result$mat_A %*% estimated_weights) / length(A)  - result$B)
    #print(estimated_weights)
    return(estimated_weights)
  }
  
  
}




neg_objective_func<-function(lam, mat_A,B){
  return(-mean(rho(t(lam) %*% mat_A)) + lam%*%B)
}


cal_solve_weights_separate<-function(str_num, str_indicator, A, pai, hX_1, hX_0,standard_weights=TRUE,EL=FALSE){
  
  n = length(A)
  weights = rep(1,n)
  
  for(i in 1:str_num){
    indistr = which(str_indicator == i)
    hX_1 = as.matrix(hX_1)
    hX_0 = as.matrix(hX_0)
    weightsistr = cal_solve_weights(1,rep(1,length(indistr)),
                                    A[indistr],
                                    pai[indistr],
                                    hX_1[indistr,], 
                                    hX_0[indistr,],
                                    standard_weights=standard_weights,
                                    EL=EL)
    weights[indistr] = weightsistr
    
  }
  return(weights)
}



cal_solve_weights_agg<-function(A, pai, hX_1, hX_0, square=FALSE){
  result = generate_objective_func_agg(A,
                                       pai,
                                       as.matrix(hX_1), 
                                       as.matrix(hX_0))
  square = TRUE
  if (square) {
    XTX = result$mat_A %*% t(result$mat_A) / length(A)
    
    lam = ginv(XTX)%*% (rowMeans(result$mat_A) - result$B)
    #print(lam)
    estimated_weights = -t(t(lam) %*% result$mat_A) + 1
    #print(t(estimated_weights))
    #print((result$mat_A %*% estimated_weights) / length(A)  - result$B)
    return(estimated_weights)
  }
  
  
  
  
  
  
  start_values = rep(0,length(result$B))
  opti_result = optim(par = start_values, 
                      fn = function(x) neg_objective_func(x, result$mat_A,result$B))
  
  lam_optimal = opti_result$par
  #print(lam_optimal)
  estimated_weights = t(d.rho(t(lam_optimal) %*% result$mat_A))
  #print(max(estimated_weights))
  #print(min(estimated_weights))
  

  return(estimated_weights)
  
}





