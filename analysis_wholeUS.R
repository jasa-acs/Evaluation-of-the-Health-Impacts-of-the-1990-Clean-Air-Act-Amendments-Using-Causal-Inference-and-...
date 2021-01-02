set.seed(2)

nboot<-50
calp<-0.1

source('functions/exmatch_trt_v2.R')
source('functions/mdmatch_conf_v2.R')
source('functions/match_rn_app_v2.R')
source('functions/match_rn_boot6.R')

## load the data for analysis ##
load('analysis_data.RData')

## for both years of data ##
for (yr in 1:2){
  
  ## pick off that year's dataset ##
  adat<-mdat[[yr]]
  
  N<-nrow(adat)

  ## create separate matricies of data for the observed exposures, counterfactual exposures, and confounders ##
  Eobs<-as.matrix(adat[,c('pmWith','ozWith')])
  Ecf<-as.matrix(adat[,c('pmNo','ozNo')])
  X<-as.matrix(adat[,grep('X',names(adat))])
  
  P<-ncol(X)
  
  ## create an external text file to save results ##
  sink(paste0('results_yr',yr,'.txt'))
  
  cat('=====================================\n')
  cat("Data Features\n")
  cat('=====================================\n')
  cat(paste0('Units in original dataset=',N,'\n'))
  
  ##############
  ## MATCHING ##
  ##############
  
  ## omega values, tolerances for exact matching on E ##
  ecut<-apply(Ecf,2,sd)*calp
  xx<-matrix(ecut,nrow=nrow(Ecf),ncol=ncol(Ecf),byrow=T)
  
  ## caliper for mahalanobis distance matching on confounders ##
  quan<-NULL
  S_inv<-solve(cov(X))
  for (j in 1:N){
    temp<-matrix(rep(X[j,],N-1),nrow=N-1,byrow=T)
    quan<-c(quan,quantile(sqrt(rowSums((temp-X[-j,])%*%S_inv*(temp-X[-j,]))),.1))
  }
  
  ## carry out matching ##
  ## the matching function outputs a dataframe with columns
  ## trtind: the indices of the units for which we found matches for their counterfactual pollutants and their confounders
  ## ctlind: the indices of the matched units for the unit given in the trtind column in the same row
  ## note that units for which multiple matches are found will be repeated in the trtind column,
  ## with each row having a different ctlind value corresponding to a different matched unit
  orig.matchout.1<-match_rn_app_v2(trtobs=Eobs,trtcf=Ecf,confounders=X,trtdiff=xx,mdqtl=mean(quan))
  orig.matchind.1<-as.data.frame(orig.matchout.1)
  names(orig.matchind.1)<-c('trtind','ctlind')
  keep.1<-unique(orig.matchind.1$trtind)
  if (yr==1){
    save(orig.matchind.1,keep.1,file='match_output.RData')
  }
  
  cat(paste0('Units in trimmed dataset=',length(keep.1),'\n'))
  
  cat('\n')
  cat('\n')
  cat('\n')
  
  cat('=====================================\n')
  cat("Matching Point Estimates\n")
  cat('=====================================\n')
  
  pname<-c('Mortality Events','Neuro Events','CVD Events')
  
  ## bias correction and point estimates ##
  allest<-NULL
  
  fmla<-as.formula(paste0('y~ s(pm25) + s(o3) + ',paste0('s(conf',1:P,')',collapse='+')))
  fmla_lm<-as.formula(paste0('y~ pm25 + o3 + ',paste0('conf',1:P,collapse='+')))
  
  ## for each of the 3 health outcomes ##
  for (i in 2:4){
    ## pick off the health outcome ##
    Yobs<-as.matrix(adat[,i])
    
    ## fit regression for bias correction ##
    regdat<-data.frame(Yobs,Eobs,X)
    names(regdat)<-c('y','pm25','o3',paste0('conf',1:P))
    fit.1<-bam(fmla, family="poisson", data=regdat,discrete=T)
    ## for each CF level in matched data, use regression function to predict outcome with true confounders and with confounders of matches ##
    predXi<-data.frame(Ecf[orig.matchind.1$trtind,],X[orig.matchind.1$trtind,])
    predXj<-data.frame(Ecf[orig.matchind.1$trtind,],X[orig.matchind.1$ctlind,])
    names(predXi)<-c('pm25','o3',paste0('conf',1:P))
    names(predXj)<-c('pm25','o3',paste0('conf',1:P))
    ## take difference to be used for bias correction ##
    diffmeans<-predict(fit.1, predXi, type="response")-predict(fit.1, predXj, type="response")
    ## impute counterfactuals ##
    EYcf<-tapply(Yobs[orig.matchind.1$ctlind]+diffmeans,orig.matchind.1$trtind,mean)
    ## difference in imputed counterfactuals and observed number of deaths in each area ## 
    orig.est.1<-sum(EYcf-Yobs[as.numeric(names(EYcf))])
    
    ## write results to the text file ##
    cat(paste0('Number of ',pname[i-1]," Prevented: ",round(orig.est.1,4),'\n'))
    allest<-c(allest,orig.est.1)
  }
  
  
  cat('\n')
  cat('\n')
  cat('\n')
  
  cat('=====================================\n')
  cat("All Matching Results\n")
  cat('=====================================\n')

  ######################
  ## do bootstrapping ##
  ######################
  
  bootind1<-list()
  bootind2<-list()
  for (j in 1:nboot){
    
    ## bootstrap from the trimmed sample ##
    bootind1<-c(bootind1,list(sample(x=keep.1,size=length(keep.1),replace=T)))
    
    ## bootstrap from the untrimmed sample ##
    bootind2<-c(bootind2,list(sample(x=1:N,size=N,replace=T)))
    
  }
  
  ## run the matching on the bootstrapped data as described in Section 3 of the paper ##
  bootmatch<-mapply(bootind1,bootind2,FUN=match_rn_boot6,MoreArgs=list(Eobs=Eobs,Ecf=Ecf,X=X,matchind=orig.matchind.1),SIMPLIFY=T)
  
  ##########################################################
  ## compute matching estimates for each bootstrap sample ##
  ##########################################################
  
  ## create an object to save results ##
  boot_results<-NULL
  
  ## for each health outcome ##
  for (i in 2:4){
    Yobs<-as.matrix(adat[,i])
    
    est.boot<-NULL
    for (j in 1:nboot){
      
      bootmatchj<-bootmatch[[j]]
      
      ## fit GAM with matched data ##
      regdat<-data.frame(Yobs[bootind2[[j]]],Eobs[bootind2[[j]],],X[bootind2[[j]],])
      names(regdat)<-c('y','pm25','o3',paste0('conf',1:P))
      fit.1<-bam(fmla, family="poisson", data=regdat,discrete=T)
      ## for each CF level in matched data, use regression function to predict outcome with true confounders and with confounders of matches ##
      predXi<-data.frame(Ecf[bootmatchj$trtind,],X[bootmatchj$trtind,])
      predXj<-data.frame(Ecf[bootmatchj$trtind,],X[bootmatchj$ctlind,])
      names(predXi)<-c('pm25','o3',paste0('conf',1:P))
      names(predXj)<-c('pm25','o3',paste0('conf',1:P))
      ## take difference to be used for bias correction ##
      diffmeans<-predict(fit.1, predXi, type="response")-predict(fit.1, predXj, type="response")
      ## impute counterfactuals ##
      EYcf<-tapply(Yobs[bootmatchj$ctlind]+diffmeans,bootmatchj$unqind,mean)
      
      ## difference in imputed counterfactuals and observed number of deaths in each area ## 
      ice_hat.1<-EYcf-Yobs[bootind1[[j]]]
      ## sum over all areas to get total number of deaths prevented ##
      est.boot<-c(est.boot,sum(ice_hat.1))
    }
    
    boot_results<-rbind(boot_results,est.boot)
  }
  

  ## write results to the text file ##
  for (i in 1:3){
    cat(paste0(pname[i],": ",round(allest[i],4)," (",round(quantile(boot_results[i,],.025),4),',',round(quantile(boot_results[i,],.975),4),')\n'))
  }
  
  cat('\n')
  cat('\n')
  cat('\n')
  
  cat('=====================================\n')
  cat("BART Results\n")
  cat('=====================================\n')
  
  ##########
  ## BART ##
  ##########
  
  ## for each health outcome ##
  for (i in 2:4){
    ## pick off the health outcome ##
    Yobs<-as.matrix(adat[,i])
    
    ## fit the BART model to the observed data ##
    ## use BART to predict counterfactuals ##
    bartfit<-bart(x.train=cbind(Eobs,X),y.train=as.numeric(Yobs),x.test=cbind(Ecf,X)[keep.1,],nskip=1000,ndpost=800,verbose=F)
    
    ## use BART to estimate the TEA (note that we are only using the units retained after trimming for estimation) ##
    bartest_tr<-sum(bartfit$yhat.test.mean-Yobs[keep.1])
    
    ## compute the credible interval ##
    tau_pd_tr<-NULL
    for (j in 1:800){
      tau_pd_tr<-c(tau_pd_tr,sum(rnorm(n=length(keep.1),mean=bartfit$yhat.test[j,],sd=bartfit$sigma[j])-Yobs[keep.1]))
    }
    bart_ci<-quantile(tau_pd_tr,c(.025,.975))
    
    ## write results to the text file ##
    cat(paste0(pname[i-1],': ',round(bartest_tr,4)," (",round(bart_ci[1],4),',',round(bart_ci[2],4),')\n'))
    
  }
  
  cat('\n')
  cat('\n')
  cat('\n')
  
  cat('=====================================\n')
  cat("Linear Model Results\n")
  cat('=====================================\n')
  
  ##################
  ## LINEAR MODEL ##
  ##################
  
  ## for each health outcome ##
  for (i in 2:4){
    ## pick off the health outcome ##
    Yobs<-as.matrix(adat[,i])
    
    ## fit linear model including all variables from unmatched data ##
    regdat<-data.frame(Yobs,Eobs,X)
    names(regdat)<-c('y','pm25','o3',paste0('conf',1:P))
    fit.3<-glm(fmla_lm, family="poisson", data=regdat)
    ## predict at counterfactual pollution levels ##
    predXi<-data.frame(Ecf[keep.1,],X[keep.1,])
    names(predXi)<-c('pm25','o3',paste0('conf',1:P))
    EYcf<-predict(fit.3, predXi, type="response")
    ice_hat.3<-EYcf-Yobs[keep.1]
    lmest<-sum(ice_hat.3)
    ## compute var-cov matrix for the predictions so that we can get the variance of the their sum ##
    vdat<-data.frame(Yobs[keep.1],Ecf[keep.1,],X[keep.1,])
    names(vdat)<-c('y','pm25','o3',paste0('conf',1:P))
    foo<-model.matrix(fmla_lm,data=vdat)
    pred_vcov<-(foo*matrix(rep(EYcf,ncol(foo)),nrow=length(keep.1),ncol=ncol(foo),byrow=F))%*%vcov(fit.3)%*%
      t(foo*matrix(rep(EYcf,ncol(foo)),nrow=length(keep.1),ncol=ncol(foo),byrow=F))
    
    ## compute the confidence limits ##
    lmlb<-sum(ice_hat.3)-1.96*sqrt(matrix(1,ncol=length(keep.1),nrow=1)%*%pred_vcov%*%matrix(1,ncol=1,nrow=length(keep.1)))
    lmub<-sum(ice_hat.3)+1.96*sqrt(matrix(1,ncol=length(keep.1),nrow=1)%*%pred_vcov%*%matrix(1,ncol=1,nrow=length(keep.1)))
    
    ## write results to the text file ##
    cat(paste0(pname[i-1],': ',round(lmest,4)," (",round(lmlb,4),',',round(lmub,4),')\n'))
  }
  
  ## stop writing to the text file ##
  sink()
  
}