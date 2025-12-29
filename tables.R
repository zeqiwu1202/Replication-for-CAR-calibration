library(calculus)
trueval<-function(Model,alpha0,alpha1,betavec0,betavec1){
  if(Model == "Model1"||Model == "Model5"){
    tv = alpha1-alpha0+(betavec1[1]-betavec0[1])*3/7+(betavec1[4]-betavec0[4])*3.8
  }
  else if(Model == "Model2"||Model == "Model6"){
    m201e<-function(x){
      60*(betavec0[1]*log(x+1)+betavec0[2]*x^2)*x^2*(1-x)^3
    }
    m202e<-function(x){
      (betavec0[3]*exp(x)+betavec0[4]/(x+3))/4
    }
    m211e<-function(x){
      60*x^2*(1-x)^3*(betavec1[1]*exp(x+2)+betavec1[2]/(x+1))
    }
    m212e<-function(x){
      betavec1[3]*x^2/4
    }
    tv = alpha1 - alpha0 + integrate(m211e,lower = 0,upper = 1)$value + 
      integrate(m212e,lower = -2,upper = 2)$value - integrate(m201e,lower = 0,upper = 1)$value -
      integrate(m202e,lower = -2,upper = 2)$value
  }
  else if(Model == "Model3"||Model == "Model7"){
    m301e<-function(x,y){
      x*y/(x+y+2)*dbeta(x,3,4)*dunif(y,-2,2)
    }
    m302e<-function(x,y,z){
      x^2*(y+z)*dbeta(x,3,4)*dunif(y,-2,2)*dnorm(z,0,1)
    }
    m311e<-function(x,y){
      (x+y)*dunif(x,-2,2)*dunif(y,0,2)
    }
    m312e<-function(x,y){
      x^2/(exp(y+2))*dunif(x,-2,2)*dbeta(y,3,4)
    }
    tv = alpha1 - alpha0 + betavec1[1]*integral(m311e,bounds = list(x=c(-2,2),y=c(0,2)))$value+
      betavec1[2]*integral(m312e,bounds = list(x=c(-2,2),y=c(0,1)))$value -
      betavec0[1]*integral(m301e,bounds = list(x=c(0,1),y=c(-2,2)))$value -
      betavec0[2]*integral(m302e,bounds = list(x=c(0,1),y=c(-2,2),z=c(-Inf,Inf)))$value
  }
  else if(Model == "Model4"||Model == "Model8"){
    m40e<-function(x){
      60*betavec0[3]*log(x+1)*x^2*(1-x)^3
    }
    m41e<-function(x){
      betavec1[3]*exp(x)/4
    }
    tv = (integrate(m41e,lower = -2,upper = 2)$value - integrate(m40e,lower = 0,upper = 1)$value)*0.5 + alpha1 - alpha0
  }
  return(tv)
}

simtable.one<-function(n,alpha0,alpha1,betavec0,betavec1,pi,Model,Iternum,method,...){
  rfre = matrix(0,nrow = 2, ncol = Iternum)
  rfcp = rep(0,Iternum)
  nsre = matrix(0,nrow = 2, ncol = Iternum)
  nscp = rep(0,Iternum)
  ns2re = matrix(0,nrow = 2, ncol = Iternum)
  ns2cp = rep(0,Iternum)
  linre = matrix(0,nrow = 2, ncol = Iternum)
  lincp = rep(0,Iternum)
  kernelre = matrix(0,nrow = 2, ncol = Iternum)
  kernelcp = rep(0,Iternum)
  gbmre = matrix(0,nrow = 2, ncol = Iternum)
  gbmcp = rep(0,Iternum)
  glmnetre = matrix(0,nrow = 2, ncol = Iternum)
  glmnetcp = rep(0,Iternum)
  pcaNNetre = matrix(0,nrow = 2, ncol = Iternum)
  pcaNNetcp = rep(0,Iternum)
  rfre_strat = matrix(0,nrow = 2, ncol = Iternum)
  rfcp_strat = rep(0,Iternum)
  nsre_strat = matrix(0,nrow = 2, ncol = Iternum)
  nscp_strat = rep(0,Iternum)
  ns2re_strat = matrix(0,nrow = 2, ncol = Iternum)
  ns2cp_strat = rep(0,Iternum)
  linre_strat = matrix(0,nrow = 2, ncol = Iternum)
  lincp_strat = rep(0,Iternum)
  kernelre_strat = matrix(0,nrow = 2, ncol = Iternum)
  kernelcp_strat = rep(0,Iternum)
  gbmre_strat = matrix(0,nrow = 2, ncol = Iternum)
  gbmcp_strat = rep(0,Iternum)
  glmnetre_strat = matrix(0,nrow = 2, ncol = Iternum)
  glmnetcp_strat = rep(0,Iternum)
  pcaNNetre_strat = matrix(0,nrow = 2, ncol = Iternum)
  pcaNNetcp_strat = rep(0,Iternum)
  FUN = match.fun(paste(paste(Model,"_",sep = ''),method,sep = ''))
  tv = trueval(Model,alpha0,alpha1,betavec0,betavec1)
  if(Model == "Model5"||Model == "Model6"||Model == "Model7"||Model == "Model8"){
    rfre_high = matrix(0,nrow = 2, ncol = Iternum)
    rfcp_high = rep(0,Iternum)
    gbmre_high = matrix(0,nrow = 2, ncol = Iternum)
    gbmcp_high = rep(0,Iternum)
    pcaNNetre_high = matrix(0,nrow = 2, ncol = Iternum)
    pcaNNetcp_high = rep(0,Iternum)
    rfre_strat_high = matrix(0,nrow = 2, ncol = Iternum)
    rfcp_strat_high = rep(0,Iternum)
    gbmre_strat_high = matrix(0,nrow = 2, ncol = Iternum)
    gbmcp_strat_high = rep(0,Iternum)
    pcaNNetre_strat_high = matrix(0,nrow = 2, ncol = Iternum)
    pcaNNetcp_strat_high = rep(0,Iternum)
    for(i in 1:Iternum){
      data = FUN(n,alpha0,alpha1,betavec0,betavec1,pi,...)
      
      rfre_temp = est_rf(data)
      nsre_temp = est_ns(data,Model)
      linre_temp = est_lin(data)
      gbmre_temp = est_gbm(data)
      glmnetre_temp = est_glmnet(data)
      #pcaNNetre_temp = est_pcaNNet(data)
      rfre_temp_high = est_rf_high(data)
      gbmre_temp_high = est_gbm_high(data)
      #pcaNNetre_temp_high = est_pcaNNet_high(data)
      
      rfre[,i] = rfre_temp[1,]
      nsre[,i] = nsre_temp[1,]
      linre[,i] = linre_temp[1,]
      gbmre[,i] = gbmre_temp[1,]
      glmnetre[,i] = glmnetre_temp[1,]
      #pcaNNetre[,i] = pcaNNetre_temp[1,]
      gbmre_high[,i] = gbmre_temp_high[1,]
      #pcaNNetre_high[,i] = pcaNNetre_temp_high[1,]
      
      rfre_strat[,i] = rfre_temp[2,]
      nsre_strat[,i] = nsre_temp[2,]
      linre_strat[,i] = linre_temp[2,]
      gbmre_strat[,i] = gbmre_temp[2,]
      glmnetre_strat[,i] = glmnetre_temp[2,]
      #pcaNNetre_strat[,i] = pcaNNetre_temp[2,]
      gbmre_strat_high[,i] = gbmre_temp_high[2,]
      #pcaNNetre_strat_high[,i] = pcaNNetre_temp_high[2,]
      
      
      if((tv >= (rfre[1,i]-1.96*rfre[2,i])) & (tv < (rfre[1,i]+1.96*rfre[2,i]))){
        rfcp[i] = 1
      }
      if((tv >= (nsre[1,i]-1.96*nsre[2,i])) & (tv < (nsre[1,i]+1.96*nsre[2,i]))){
        nscp[i] = 1
      }
      if((tv >= (linre[1,i]-1.96*linre[2,i])) & (tv < (linre[1,i]+1.96*linre[2,i]))){
        lincp[i] = 1
      }
      if((tv >= (gbmre[1,i]-1.96*gbmre[2,i])) & (tv < (gbmre[1,i]+1.96*gbmre[2,i]))){
        gbmcp[i] = 1
      }
      if((tv >= (glmnetre[1,i]-1.96*glmnetre[2,i])) & (tv < (glmnetre[1,i]+1.96*glmnetre[2,i]))){
        glmnetcp[i] = 1
      }
      #if((tv >= (pcaNNetre[1,i]-1.96*pcaNNetre[2,i])) & (tv < (pcaNNetre[1,i]+1.96*pcaNNetre[2,i]))){
      #  pcaNNetcp[i] = 1
      #}
      if((tv >= (gbmre_high[1,i]-1.96*gbmre_high[2,i])) & (tv < (gbmre_high[1,i]+1.96*gbmre_high[2,i]))){
        gbmcp_high[i] = 1
      }
      if((tv >= (rfre_high[1,i]-1.96*rfre_high[2,i])) & (tv < (rfre_high[1,i]+1.96*rfre_high[2,i]))){
        rfcp_high[i] = 1
      }
      #if((tv >= (pcaNNetre_high[1,i]-1.96*pcaNNetre_high[2,i])) & (tv < (pcaNNetre_high[1,i]+1.96*pcaNNetre_high[2,i]))){
      #  pcaNNetcp_high[i] = 1
      #}
      
      
      if((tv >= (rfre_strat[1,i]-1.96*rfre_strat[2,i])) & (tv < (rfre_strat[1,i]+1.96*rfre_strat[2,i]))){
        rfcp_strat[i] = 1
      }
      if((tv >= (nsre_strat[1,i]-1.96*nsre_strat[2,i])) & (tv < (nsre_strat[1,i]+1.96*nsre_strat[2,i]))){
        nscp_strat[i] = 1
      }
      if((tv >= (linre_strat[1,i]-1.96*linre_strat[2,i])) & (tv < (linre_strat[1,i]+1.96*linre_strat[2,i]))){
        lincp_strat[i] = 1
      }
      if((tv >= (gbmre_strat[1,i]-1.96*gbmre_strat[2,i])) & (tv < (gbmre_strat[1,i]+1.96*gbmre_strat[2,i]))){
        gbmcp_strat[i] = 1
      }
      if((tv >= (glmnetre_strat[1,i]-1.96*glmnetre_strat[2,i])) & (tv < (glmnetre_strat[1,i]+1.96*glmnetre_strat[2,i]))){
        glmnetcp_strat[i] = 1
      }
      #if((tv >= (pcaNNetre_strat[1,i]-1.96*pcaNNetre_strat[2,i])) & (tv < (pcaNNetre_strat[1,i]+1.96*pcaNNetre_strat[2,i]))){
      #  pcaNNetcp_strat[i] = 1
      #}
      if((tv >= (gbmre_strat_high[1,i]-1.96*gbmre_strat_high[2,i])) & (tv < (gbmre_strat_high[1,i]+1.96*gbmre_strat_high[2,i]))){
        gbmcp_strat_high[i] = 1
      }
      if((tv >= (rfre_strat[1,i]-1.96*rfre_strat_high[2,i])) & (tv < (rfre_strat_high[1,i]+1.96*rfre_strat_high[2,i]))){
        rfcp_strat_high[i] = 1
      }
      #if((tv >= (pcaNNetre_strat[1,i]-1.96*pcaNNetre_strat_high[2,i])) & (tv < (pcaNNetre_strat_high[1,i]+1.96*pcaNNetre_strat_high[2,i]))){
      #  pcaNNetcp_strat_high[i] = 1
      #}
      
    }
    comtabs = matrix(0,9,4)
    comcps = rep(0,9)
    comtabs[1,] = c(mean(rfre[1,]) - tv,sd(rfre[1,]),mean(rfre[2,]),mean(rfre[3,]))
    comcps[1] = sum(rfcp)/Iternum
    comtabs[2,] = c(mean(rfre_high[1,]) - tv,sd(rfre_high[1,]),mean(rfre_high[2,]),mean(rfre_high[3,]))
    comcps[2] = sum(rfcp_high)/Iternum
    comtabs[3,] = c(mean(nsre[1,]) - tv,sd(nsre[1,]),mean(nsre[2,]),mean(nsre[3,]))
    comcps[3] = sum(nscp)/Iternum
    comtabs[4,] = c(mean(linre[1,]) - tv,sd(linre[1,]),mean(linre[2,]),mean(linre[3,]))
    comcps[4] = sum(lincp)/Iternum
    comtabs[5,] = c(mean(glmnetre[1,]) - tv,sd(glmnetre[1,]),mean(glmnetre[2,]),mean(glmnetre[3,]))
    comcps[5] = sum(glmnetcp)/Iternum
    comtabs[6,] = c(mean(gbmre[1,]) - tv,sd(gbmre[1,]),mean(gbmre[2,]),mean(gbmre[3,]))
    comcps[6] = sum(gbmcp)/Iternum
    comtabs[7,] = c(mean(gbmre_high[1,]) - tv,sd(gbmre_high[1,]),mean(gbmre_high[2,]),mean(gbmre_high[3,]))
    comcps[7] = sum(gbmcp_high)/Iternum
    #comtabs[8,] = c(mean(pcaNNetre[1,]) - tv,sd(pcaNNetre[1,]),mean(pcaNNetre[2,]),mean(pcaNNetre[3,]))
    #comcps[8] = sum(pcaNNetcp)/Iternum
    #comtabs[9,] = c(mean(pcaNNetre_high[1,]) - tv,sd(pcaNNetre_high[1,]),mean(pcaNNetre_high[2,]),mean(pcaNNetre_high[3,]))
    #comcps[9] = sum(pcaNNetcp_high)/Iternum
    rownames(comtabs) = c("randforest","randomforest_high","naturalspline","linear","Elastic Net",
                          "Boosted Trees","Boosted Trees_high","pcaNNet","pcaNNet_high")
    colnames(comtabs) = c("bias","sd","se","hr")
    names(comcps) = c("randforest","randomforest_high","naturalspline","linear","Elastic Net",
                      "Boosted Trees","Boosted Trees_high","pcaNNet","pcaNNet_high")
    
    comtabs_strat = matrix(0,9,4)
    comcps_strat = rep(0,9)
    comtabs_strat[1,] = c(mean(rfre_strat[1,]) - tv,sd(rfre_strat[1,]),mean(rfre_strat[2,]),mean(rfre_strat[3,]))
    comcps_strat[1] = sum(rfcp_strat)/Iternum
    comtabs_strat[2,] = c(mean(rfre_strat_high[1,]) - tv,sd(rfre_strat_high[1,]),mean(rfre_strat_high[2,]),mean(rfre_strat_high[3,]))
    comcps_strat[2] = sum(rfcp_strat_high)/Iternum
    comtabs_strat[3,] = c(mean(nsre_strat[1,]) - tv,sd(nsre_strat[1,]),mean(nsre_strat[2,]),mean(nsre_strat[3,]))
    comcps_strat[3] = sum(nscp_strat)/Iternum
    comtabs_strat[4,] = c(mean(linre_strat[1,]) - tv,sd(linre_strat[1,]),mean(linre_strat[2,]),mean(linre_strat[3,]))
    comcps_strat[4] = sum(lincp_strat)/Iternum
    comtabs_strat[5,] = c(mean(glmnetre_strat[1,]) - tv,sd(glmnetre_strat[1,]),mean(glmnetre_strat[2,]),mean(glmnetre_strat[3,]))
    comcps_strat[5] = sum(glmnetcp_strat)/Iternum
    comtabs_strat[6,] = c(mean(gbmre_strat[1,]) - tv,sd(gbmre_strat[1,]),mean(gbmre_strat[2,]),mean(gbmre_strat[3,]))
    comcps_strat[6] = sum(gbmcp_strat)/Iternum
    comtabs_strat[7,] = c(mean(gbmre_strat_high[1,]) - tv,sd(gbmre_strat_high[1,]),mean(gbmre_strat_high[2,]),mean(gbmre_strat_high[3,]))
    comcps_strat[7] = sum(gbmcp_strat_high)/Iternum
    #comtabs_strat[8,] = c(mean(pcaNNetre_strat[1,]) - tv,sd(pcaNNetre_strat[1,]),mean(pcaNNetre_strat[2,]),mean(pcaNNetre_strat[3,]))
    #comcps_strat[8] = sum(pcaNNetcp_strat)/Iternum
    #comtabs_strat[9,] = c(mean(pcaNNetre_strat_high[1,]) - tv,sd(pcaNNetre_strat_high[1,]),mean(pcaNNetre_strat_high[2,]),mean(pcaNNetre_strat_high[3,]))
    #comcps_strat[9] = sum(pcaNNetcp_strat_high)/Iternum
    rownames(comtabs_strat) = c("randforest","randomforest_high","naturalspline","linear","Elastic Net",
                                "Boosted Trees","Boosted Trees_high","pcaNNet","pcaNNet_high")
    colnames(comtabs_strat) = c("bias","sd","se","hr")
    names(comcps_strat) = c("randforest","randomforest_high","naturalspline","linear","Elastic Net",
                            "Boosted Trees","Boosted Trees_high","pcaNNet","pcaNNet_high")
    }
  else{
    for(i in 1:Iternum){
      data = FUN(n,alpha0,alpha1,betavec0,betavec1,pi,...)
      
      rfre_temp = est_rf(data)
      nsre_temp = est_ns(data,Model)
      #ns2re_temp = est_ns2(data,Model)
      linre_temp = est_lin(data)
      kernelre_temp = est_kernel(data)
      #gbmre_temp = est_gbm(data)
      #glmnetre_temp = est_glmnet(data)
      #pcaNNetre_temp = est_pcaNNet(data)
      
      
      rfre[,i] = rfre_temp[1,]
      nsre[,i] = nsre_temp[1,]
      #ns2re[,i] = ns2re_temp[1,]
      linre[,i] = linre_temp[1,]
      kernelre[,i] = kernelre_temp[1,]
      #gbmre[,i] = gbmre_temp[1,]
      #glmnetre[,i] = glmnetre_temp[1,]
      #pcaNNetre[,i] = pcaNNetre_temp[1,]
      
      rfre_strat[,i] = rfre_temp[2,]
      nsre_strat[,i] = nsre_temp[2,]
      #ns2re_strat[,i] = ns2re_temp[2,]
      linre_strat[,i] = linre_temp[2,]
      kernelre_strat[,i] = kernelre_temp[2,]
      #gbmre_strat[,i] = gbmre_temp[2,]
      #glmnetre_strat[,i] = glmnetre_temp[2,]
      #pcaNNetre_strat[,i] = pcaNNetre_temp[2,]
      
      
      if((tv >= (rfre[1,i]-1.96*rfre[2,i])) & (tv < (rfre[1,i]+1.96*rfre[2,i]))){
        rfcp[i] = 1
      }
      if((tv >= (nsre[1,i]-1.96*nsre[2,i])) & (tv < (nsre[1,i]+1.96*nsre[2,i]))){
        nscp[i] = 1
      }
      #if((tv >= (ns2re[1,i]-1.96*ns2re[2,i])) & (tv < (ns2re[1,i]+1.96*ns2re[2,i]))){
      #  ns2cp[i] = 1
      #}
      if((tv >= (linre[1,i]-1.96*linre[2,i])) & (tv < (linre[1,i]+1.96*linre[2,i]))){
        lincp[i] = 1
      }
      if((tv >= (kernelre[1,i]-1.96*kernelre[2,i])) & (tv < (kernelre[1,i]+1.96*kernelre[2,i]))){
        kernelcp[i] = 1
      }
      #if((tv >= (gbmre[1,i]-1.96*gbmre[2,i])) & (tv < (gbmre[1,i]+1.96*gbmre[2,i]))){
      #  gbmcp[i] = 1
      #}
      #if((tv >= (glmnetre[1,i]-1.96*glmnetre[2,i])) & (tv < (glmnetre[1,i]+1.96*glmnetre[2,i]))){
      #  glmnetcp[i] = 1
      #}
      #if((tv >= (pcaNNetre[1,i]-1.96*pcaNNetre[2,i])) & (tv < (pcaNNetre[1,i]+1.96*pcaNNetre[2,i]))){
      #  pcaNNetcp[i] = 1
      #}
      
      
      if((tv >= (rfre_strat[1,i]-1.96*rfre_strat[2,i])) & (tv < (rfre_strat[1,i]+1.96*rfre_strat[2,i]))){
        rfcp_strat[i] = 1
      }
      if((tv >= (nsre_strat[1,i]-1.96*nsre_strat[2,i])) & (tv < (nsre_strat[1,i]+1.96*nsre_strat[2,i]))){
        nscp_strat[i] = 1
      }
      #if((tv >= (ns2re_strat[1,i]-1.96*ns2re_strat[2,i])) & (tv < (ns2re_strat[1,i]+1.96*ns2re_strat[2,i]))){
      #  ns2cp_strat[i] = 1
      #}
      if((tv >= (linre_strat[1,i]-1.96*linre_strat[2,i])) & (tv < (linre_strat[1,i]+1.96*linre_strat[2,i]))){
        lincp_strat[i] = 1
      }
      if((tv >= (kernelre_strat[1,i]-1.96*kernelre_strat[2,i])) & (tv < (kernelre_strat[1,i]+1.96*kernelre_strat[2,i]))){
        kernelcp_strat[i] = 1
      }
      #if((tv >= (gbmre_strat[1,i]-1.96*gbmre_strat[2,i])) & (tv < (gbmre_strat[1,i]+1.96*gbmre_strat[2,i]))){
      #  gbmcp_strat[i] = 1
      #}
      #if((tv >= (glmnetre_strat[1,i]-1.96*glmnetre_strat[2,i])) & (tv < (glmnetre_strat[1,i]+1.96*glmnetre_strat[2,i]))){
      #  glmnetcp_strat[i] = 1
      #}
      #if((tv >= (pcaNNetre_strat[1,i]-1.96*pcaNNetre_strat[2,i])) & (tv < (pcaNNetre_strat[1,i]+1.96*pcaNNetre_strat[2,i]))){
      #  pcaNNetcp_strat[i] = 1
      #}
    }
    comtabs = matrix(0,8,3)
    comcps = rep(0,8)
    comtabs[1,] = c(mean(rfre[1,]) - tv,sd(rfre[1,]),mean(rfre[2,]))
    comcps[1] = sum(rfcp)/Iternum
    comtabs[2,] = c(mean(nsre[1,]) - tv,sd(nsre[1,]),mean(nsre[2,]))
    comcps[2] = sum(nscp)/Iternum
    #comtabs[3,] = c(mean(ns2re[1,]) - tv,sd(ns2re[1,]),mean(ns2re[2,]))
    #comcps[3] = sum(ns2cp)/Iternum
    comtabs[4,] = c(mean(linre[1,]) - tv,sd(linre[1,]),mean(linre[2,]))
    comcps[4] = sum(lincp)/Iternum
    comtabs[5,] = c(mean(kernelre[1,]) - tv,sd(kernelre[1,]),mean(kernelre[2,]))
    comcps[5] = sum(kernelcp)/Iternum
    #comtabs[6,] = c(mean(glmnetre[1,]) - tv,sd(glmnetre[1,]),mean(glmnetre[2,]))
    #comcps[6] = sum(glmnetcp)/Iternum
    #comtabs[7,] = c(mean(gbmre[1,]) - tv,sd(gbmre[1,]),mean(gbmre[2,]))
    #comcps[7] = sum(gbmcp)/Iternum
    #comtabs[8,] = c(mean(pcaNNetre[1,]) - tv,sd(pcaNNetre[1,]),mean(pcaNNetre[2,]))
    #comcps[8] = sum(pcaNNetcp)/Iternum
    rownames(comtabs) = c("randforest","naturalspline","natual spline with interaction","linear","kernel","Elastic Net","Boosted Trees","pcaNNet")
    colnames(comtabs) = c("bias","sd","se")
    names(comcps) = c("randforest","naturalspline","natual spline with interaction","linear","kernel","Elastic Net","Boosted Trees","pcaNNet")
    
    comtabs_strat = matrix(0,8,3)
    comcps_strat = rep(0,8)
    comtabs_strat[1,] = c(mean(rfre_strat[1,]) - tv,sd(rfre_strat[1,]),mean(rfre_strat[2,]))
    comcps_strat[1] = sum(rfcp_strat)/Iternum
    comtabs_strat[2,] = c(mean(nsre_strat[1,]) - tv,sd(nsre_strat[1,]),mean(nsre_strat[2,]))
    comcps_strat[2] = sum(nscp_strat)/Iternum
    #comtabs_strat[3,] = c(mean(ns2re_strat[1,]) - tv,sd(ns2re_strat[1,]),mean(ns2re_strat[2,]))
    #comcps_strat[3] = sum(ns2cp_strat)/Iternum
    comtabs_strat[4,] = c(mean(linre_strat[1,]) - tv,sd(linre_strat[1,]),mean(linre_strat[2,]))
    comcps_strat[4] = sum(lincp_strat)/Iternum
    comtabs_strat[5,] = c(mean(kernelre_strat[1,]) - tv,sd(kernelre_strat[1,]),mean(kernelre_strat[2,]))
    comcps_strat[5] = sum(kernelcp_strat)/Iternum
    #comtabs_strat[6,] = c(mean(glmnetre_strat[1,]) - tv,sd(glmnetre_strat[1,]),mean(glmnetre_strat[2,]))
    #comcps_strat[6] = sum(glmnetcp_strat)/Iternum
    #comtabs_strat[7,] = c(mean(gbmre_strat[1,]) - tv,sd(gbmre_strat[1,]),mean(gbmre_strat[2,]))
    #comcps_strat[7] = sum(gbmcp_strat)/Iternum
    #comtabs_strat[8,] = c(mean(pcaNNetre_strat[1,]) - tv,sd(pcaNNetre_strat[1,]),mean(pcaNNetre_strat[2,]))
    #comcps_strat[8] = sum(pcaNNetcp_strat)/Iternum
    rownames(comtabs_strat) = c("randforest","natural spline","natual spline with interaction"," linear","kernel","Elastic Net","Boosted Trees","pcaNNet")
    colnames(comtabs_strat) = c("bias","sd","se")
    names(comcps_strat) = c("randforest","naturalspline","natual spline with interaction","linear","kernel","Elastic Net","Boosted Trees","pcaNNet")
    
    }
  return(list(comtabs = comtabs, comcps = comcps, comtabs_strat = comtabs_strat, comcps_strat = comcps_strat))
}

sim_total<-function(n,Iternum){
  p = 200
  lambda = 0.75
  weight = rep(0.5,2)
  omega = c(0.1,0.1,0.4,0.4)
  sigma1 = 3
  sigma0 = 1
  
  #Model 1 
  alpha0 = 1
  alpha1 = 4

  betavec0 = c(75,35,125,80)
  betavec1 = c(100,80,60,40)
  
  #equal allocation
  pi = 1/2
  #SRSsim1 = simtable.one(n,alpha0,alpha1,betavec0,betavec1,pi,"Model1",Iternum,"SRS",sigma0,sigma1)
  #SRSsim1_5000 = simtable.one(5000,alpha0,alpha1,betavec0,betavec1,pi,"Model1",Iternum,"SRS",sigma0,sigma1)
  
  #WEIsim1 = simtable.one(n,alpha0,alpha1,betavec0,betavec1,pi,"Model1",Iternum,"WEI",sigma0,sigma1)
  
  SBRsim1 = simtable.one(n,alpha0,alpha1,betavec0,betavec1,pi,"Model1",Iternum,"SBR",sigma0,sigma1)
  
  #BCDsim1 = simtable.one(n,alpha0,alpha1,betavec0,betavec1,pi,"Model1",Iternum,"BCD",sigma0,sigma1,lambda)
  
  #PSsim1 = simtable.one(n,alpha0,alpha1,betavec0,betavec1,pi,"Model1",Iternum,"PS",sigma0,sigma1,weight,lambda)
  
  #HHsim1 = simtable.one(n,alpha0,alpha1,betavec0,betavec1,pi,"Model1",Iternum,"HH",sigma0,sigma1,omega,lambda)
  
  #unequal allocation
  #pi = 2/3
  
  #SRSsim123 = simtable.one(n,alpha0,alpha1,betavec0,betavec1,pi,"Model1",Iternum,"SRS",sigma0,sigma1)
  
  #  WEIsim123 = simtable.one(n,alpha0,alpha1,betavec0,betavec1,pi,"Model1",Iternum,"WEI",sigma0,sigma1)
  
  #SBRsim123 = simtable.one(n,alpha0,alpha1,betavec0,betavec1,pi,"Model1",Iternum,"SBR",sigma0,sigma1)
  
  #  BCDsim123 = simtable.one(n,alpha0,alpha1,betavec0,betavec1,pi,"Model1",Iternum,"BCD",sigma0,sigma1,lambda)
  
  #PSsim123 = simtable.one(n,alpha0,alpha1,betavec0,betavec1,pi,"Model1",Iternum,"PS",sigma0,sigma1,weight,lambda)
  
  #  HHsim123 = simtable.one(n,alpha0,alpha1,betavec0,betavec1,pi,"Model1",Iternum,"HH",sigma0,sigma1,omega,lambda)
  
  #pi = 4/5
  
  #SRSsim105_5000 = simtable.one(5000,alpha0,alpha1,betavec0,betavec1,pi,"Model1",Iternum,"SRS",sigma0,sigma1)
  
  #  WEIsim123 = simtable.one(n,alpha0,alpha1,betavec0,betavec1,pi,"Model1",Iternum,"WEI",sigma0,sigma1)
  
  #SBRsim105 = simtable.one(n,alpha0,alpha1,betavec0,betavec1,pi,"Model1",Iternum,"SBR",sigma0,sigma1)
  
  #  BCDsim123 = simtable.one(n,alpha0,alpha1,betavec0,betavec1,pi,"Model1",Iternum,"BCD",sigma0,sigma1,lambda)

  #Model 2
  alpha0 = -3
  alpha1 = 0
  betavec0 = c(10,24,15,20)
  betavec1 = c(20,27,10)
  #pi = 1/2
  #pi = 1/2
  #SRSsim2 = simtable.one(n,alpha0,alpha1,betavec0,betavec1,pi,"Model2",Iternum,"SRS",sigma0,sigma1)
  
  SBRsim2 = simtable.one(n,alpha0,alpha1,betavec0,betavec1,pi,"Model2",Iternum,"SBR",sigma0,sigma1)
  
  #PSsim2 = simtable.one(n,alpha0,alpha1,betavec0,betavec1,pi,"Model2",Iternum,"PS",sigma0,sigma1,weight,lambda)
  
  #pi = 2/3
  #SRSsim223 = simtable.one(n,alpha0,alpha1,betavec0,betavec1,pi,"Model2",Iternum,"SRS",sigma0,sigma1)
  
  #SBRsim223 = simtable.one(n,alpha0,alpha1,betavec0,betavec1,pi,"Model2",Iternum,"SBR",sigma0,sigma1)
  
  #PSsim223 = simtable.one(n,alpha0,alpha1,betavec0,betavec1,pi,"Model2",Iternum,"PS",sigma0,sigma1,weight,lambda)
  
  #Model 3
  alpha0 = 5
  alpha1 = 2
  betavec0 = c(24,43)
  betavec1 = c(22,35)
  
  #pi = 1/2
  #pi = 1/2
  #SRSsim3 = simtable.one(n,alpha0,alpha1,betavec0,betavec1,pi,"Model3",Iternum,"SRS",sigma0,sigma1)
  
  SBRsim3 = simtable.one(n,alpha0,alpha1,betavec0,betavec1,pi,"Model3",Iternum,"SBR",sigma0,sigma1)
  
  #PSsim3 = simtable.one(n,alpha0,alpha1,betavec0,betavec1,pi,"Model3",Iternum,"PS",sigma0,sigma1,weight,lambda)
  
  #pi = 2/3
  #SRSsim323 = simtable.one(n,alpha0,alpha1,betavec0,betavec1,pi,"Model3",Iternum,"SRS",sigma0,sigma1)
  
  #SBRsim323 = simtable.one(n,alpha0,alpha1,betavec0,betavec1,pi,"Model3",Iternum,"SBR",sigma0,sigma1)
  
  #PSsim323 = simtable.one(n,alpha0,alpha1,betavec0,betavec1,pi,"Model3",Iternum,"PS",sigma0,sigma1,weight,lambda)
  
  #Model 4
  alpha0 = 5
  alpha1 = 5
  betavec0 = c(4,6,2)
  betavec1 = c(4,6,5)
  #pi = 1/2
  #pi = 1/2
  #SRSsim4 = simtable.one(n,alpha0,alpha1,betavec0,betavec1,pi,"Model4",Iternum,"SRS",sigma0,sigma1)
  #SRSsim4_5000 = simtable.one(5000,alpha0,alpha1,betavec0,betavec1,pi,"Model4",Iternum,"SRS",sigma0,sigma1)
  SBRsim4 = simtable.one(n,alpha0,alpha1,betavec0,betavec1,pi,"Model4",Iternum,"SBR",sigma0,sigma1)
  
  #PSsim4 = simtable.one(n,alpha0,alpha1,betavec0,betavec1,pi,"Model4",Iternum,"PS",sigma0,sigma1,weight,lambda)
  
  #pi = 2/3
  #SRSsim423 = simtable.one(n,alpha0,alpha1,betavec0,betavec1,pi,"Model4",Iternum,"SRS",sigma0,sigma1)
  
  #SBRsim423 = simtable.one(n,alpha0,alpha1,betavec0,betavec1,pi,"Model4",Iternum,"SBR",sigma0,sigma1)
  
  #PSsim423 = simtable.one(n,alpha0,alpha1,betavec0,betavec1,pi,"Model4",Iternum,"PS",sigma0,sigma1,weight,lambda)
  
  #Model 5
  #alpha0 = -3
  #alpha1 = 0
  #betavec0 = c(20,34,42,10)
  #betavec1 = c(20,17,10)

  #pi = 1/2
  #SRSsim5 = simtable.one(n,alpha0,alpha1,betavec0,betavec1,pi,"Model5",Iternum,"SRS",sigma0,sigma1)
  #SRSsim5_5000 = simtable.one(5000,alpha0,alpha1,betavec0,betavec1,pi,"Model5",Iternum,"SRS",sigma0,sigma1)
  #SBRsim5 = simtable.one(n,alpha0,alpha1,betavec0,betavec1,pi,"Model5",Iternum,"SBR",sigma0,sigma1)
  
  #PSsim5 = simtable.one(n,alpha0,alpha1,betavec0,betavec1,pi,"Model5",Iternum,"PS",sigma0,sigma1,weight,lambda)
  
  #pi = 2/3
  #SRSsim523 = simtable.one(n,alpha0,alpha1,betavec0,betavec1,pi,"Model5",Iternum,"SRS",sigma0,sigma1)
  
  #SBRsim523 = simtable.one(n,alpha0,alpha1,betavec0,betavec1,pi,"Model5",Iternum,"SBR",sigma0,sigma1)
  
  #PSsim523 = simtable.one(n,alpha0,alpha1,betavec0,betavec1,pi,"Model5",Iternum,"PS",sigma0,sigma1,weight,lambda)
  
  #Model 6
  #alpha0 = -3
  #alpha1 = 0
  #betavec0 = c(10,20,15)
  #betavec1 = c(30,40)

  #pi = 1/2
  #SRSsim6 = simtable.one(n,alpha0,alpha1,betavec0,betavec1,pi,"Model6",Iternum,"SRS",sigma0,sigma1)
  
  #SBRsim6 = simtable.one(n,alpha0,alpha1,betavec0,betavec1,pi,"Model6",Iternum,"SBR",sigma0,sigma1)
  
  #PSsim6 = simtable.one(n,alpha0,alpha1,betavec0,betavec1,pi,"Model6",Iternum,"PS",sigma0,sigma1,weight,lambda)
  
  #pi = 2/3
  #SRSsim623 = simtable.one(n,alpha0,alpha1,betavec0,betavec1,pi,"Model6",Iternum,"SRS",sigma0,sigma1)
  
  #SBRsim623 = simtable.one(n,alpha0,alpha1,betavec0,betavec1,pi,"Model6",Iternum,"SBR",sigma0,sigma1)
  
  #PSsim623 = simtable.one(n,alpha0,alpha1,betavec0,betavec1,pi,"Model6",Iternum,"PS",sigma0,sigma1,weight,lambda)
  
  #Model 7
  #alpha0 = -3
  #alpha1 = 1
  #betavec0 = c(10,7,20)
  #betavec1 = c(6,20,16)
  #sigma0 = 1
  #sigma1 = 2

  #pi = 1/2
  #SRSsim7 = simtable.one(n,alpha0,alpha1,betavec0,betavec1,pi,"Model7",Iternum,"SRS",sigma0,sigma1,p)
  
  #SBRsim7 = simtable.one(n,alpha0,alpha1,betavec0,betavec1,pi,"Model7",Iternum,"SBR",sigma0,sigma1,p)
  
  #PSsim7 = simtable.one(n,alpha0,alpha1,betavec0,betavec1,pi,"Model7",Iternum,"PS",sigma0,sigma1,p,weight,lambda)
  
  #pi = 2/3
  #SRSsim723 = simtable.one(n,alpha0,alpha1,betavec0,betavec1,pi,"Model7",Iternum,"SRS",sigma0,sigma1,p)
  
  #SBRsim723 = simtable.one(n,alpha0,alpha1,betavec0,betavec1,pi,"Model7",Iternum,"SBR",sigma0,sigma1,p)
  
  #PSsim723 = simtable.one(n,alpha0,alpha1,betavec0,betavec1,pi,"Model7",Iternum,"PS",sigma0,sigma1,p,weight,lambda)
  
  #Model 8
  #alpha0 = -3
  #alpha1 = 1
  #betavec0 = c(10,7,20)
  #betavec1 = 30
  #sigma0 = 1
  #sigma1 = 2

  #pi = 1/2
  #SRSsim8 = simtable.one(n,alpha0,alpha1,betavec0,betavec1,pi,"Model8",Iternum,"SRS",sigma0,sigma1,p)
  
  #SBRsim8 = simtable.one(n,alpha0,alpha1,betavec0,betavec1,pi,"Model8",Iternum,"SBR",sigma0,sigma1,p)
  
  #PSsim8 = simtable.one(n,alpha0,alpha1,betavec0,betavec1,pi,"Model8",Iternum,"PS",sigma0,sigma1,p,weight,lambda)
  
  #pi = 2/3
  #SRSsim823 = simtable.one(n,alpha0,alpha1,betavec0,betavec1,pi,"Model8",Iternum,"SRS",sigma0,sigma1,p)
  
  #SBRsim823 = simtable.one(n,alpha0,alpha1,betavec0,betavec1,pi,"Model8",Iternum,"SBR",sigma0,sigma1,p)
  
  #PSsim823 = simtable.one(n,alpha0,alpha1,betavec0,betavec1,pi,"Model8",Iternum,"PS",sigma0,sigma1,p,weight,lambda)
  

  return(list(SBRsim1 = SBRsim1,SBRsim2 = SBRsim2,SBRsim3 = SBRsim3,SBRsim4 = SBRsim4))
}

library(xtable)
textab<-function(result){
  nums = length(result)
  tabs = NULL
  for(i in seq(1,nums,3)){
    tab_temp = NULL
    for(j in 1:3){
      l = length(result[[i+j-1]]$comcps)
      tab_temp = cbind(tab_temp,result[[i+j-1]]$comtabs[-c(3,l),],result[[i+j-1]]$comcps[-c(3,l)])
    }
    tabs = rbind(tabs,tab_temp)
  }
  return(xtable(tabs))
}

textab5000<-function(result){
  nums = length(result)
  tabs = NULL
  for(i in seq(1,nums,1)){
    tab_temp = NULL
    l = length(result[[i]]$comcps)
    tab_temp = cbind(tab_temp,result[[i]]$comtabs[-c(3,l),],result[[i]]$comcps[-c(3,l)])
    tabs = rbind(tabs,tab_temp)
  }
  return(xtable(tabs))
}

textab2<-function(result){
  nums = length(result)
  tabs = NULL
  for(i in seq(1,nums,3)){
    tab_temp = NULL
    tab_strat_temp = NULL
    tab_S_temp = NULL
    for(j in 1:3){
      tab_temp = cbind(tab_temp,result[[i+j-1]]$comtabs[c(1,2,4,5),],result[[i+j-1]]$comcps[c(1,2,4,5)])
      tab_strat_temp = cbind(tab_strat_temp,result[[i+j-1]]$comtabs_strat[c(1,2,4,5),],result[[i+j-1]]$comcps_strat[c(1,2,4,5)])
      tab_S_temp = cbind(tab_S_temp,result[[i+j-1]]$comtabs_S[c(1,2,4,5),],result[[i+j-1]]$comcps_S[c(1,2,4,5)])
    }
    for(k in 1:nrow(tab_temp))
    tabs = rbind(tabs,tab_temp[k,],tab_S_temp[k,],tab_strat_temp[k,])
  }
  return(xtable(tabs))
}

