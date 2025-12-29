createFolds_by_strat<-function(index.set, num_of_k, input_strat,input_treatment){
  strt = unique(input_strat)
  strt_num = length(strt)
  result = list()
  for (k in 1:num_of_k){
    result[paste("Fold",k,sep = "")] = c()
  }
  for(i in 1:strt_num){
    ind_strat_treat = which(input_strat == strt[i] & input_treatment == 1)
    ind_strat_control = which(input_strat == strt[i] & input_treatment == 0)
    this_strt_treat = createFolds(ind_strat_treat,num_of_k)
    this_strt_control = createFolds(ind_strat_control,num_of_k)
    for (k in 1:num_of_k){
      eval(parse(text = paste("result$Fold",k,"= c(result$Fold",k,",ind_strat_treat[this_strt_treat$Fold",k,"])",sep = "")))
      eval(parse(text = paste("result$Fold",k,"= c(result$Fold",k,",ind_strat_control[this_strt_control$Fold",k,"])",sep = "")))
    }
    

  }
  return(result)
  
  
}