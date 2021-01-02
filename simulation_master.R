## Here I have specified the simulation parameters so that
## reps=1, one dataset will be simulated, increase this to increase the number of simulated datasets
## Ngrid=70, the data are generated on a 70 x 70 cell grid so that the sample size (number of grid cells) is 70^2=4900
## nboot=50, 50 bootstrap samples will be used to produce confidence intervals for the matching estimator
## simtype=1, the data generating process for simulation 1 as described in Section 4 of the manuscript (simtype=2 and simtype=3 give the DGPs for simulations 2 and 3 in the manuscript)
## calp=0.15, setting omega=0.15 (omega defines the matching tolerances) as described in Section 4 of the manuscript (calp=0.2 and calp=0.25 are also used to produce simulations in the manuscript)

reps<-1
Ngrid<-70
nboot<-50
simtype<-1
calp<-0.15

## load required libraries ##
library(mvtnorm)
library(dplyr)
library(splines)
library(BayesTree)
library(gstat)
library(mgcv)

## load the functions needed to run the matching and bootstrap ##
source('functions/exmatch_trt_v2.R')
source('functions/mdmatch_conf_v2.R')
source('functions/match_rn_app_v2.R')
source('functions/match_rn_boot6.R')


#######################
## 1. CONSTRUCT DATA ##
#######################

## seed for data reproducibility and consistency of pollution and confounder data across simulations ##
set.seed(4)

## simulate 2 factual (observed) pollutants similar to pm2.5 annual average and o3 summer average ##
xy<-expand.grid(1:Ngrid,1:Ngrid)
names(xy)<-c('x','y')
N<-nrow(xy)

pm_dummy <- gstat(formula=z~1, locations=~x+y, dummy=T, beta=12.18, model=vgm(nugget=3,psill=5,model="Exp",range=10), nmax=20)
pm_spat <- predict(pm_dummy, newdata=xy, nsim=1)

oz_dummy <- gstat(formula=z~1, locations=~x+y, dummy=T, beta=0.05, model=vgm(nugget=.00005,psill=.00005,model="Exp",range=10), nmax=20)
oz_spat <- predict(oz_dummy, newdata=xy, nsim=1)

Eobs<-as.matrix(data.frame(pm_spat$sim1,oz_spat$sim1*1000))
colnames(Eobs)<-c('pm25','o3')

## simulate 2 counterfactual pollutants with larger mean and larger spread ##
Ecf<-Eobs+data.frame(runif(n=N,min=0,max=2.5*sd(Eobs[,1],na.rm=T)),runif(n=N,min=0,max=2.5*sd(Eobs[,2],na.rm=T)))
Ecf<-as.matrix(Ecf)
colnames(Ecf)<-c('pm25','o3')

## simulate confounders ##
betaE<-matrix(c(c(-0.2291632, -0.2212561, -0.1103875, -0.2885972,  0.1754248),c(12,-5,-3.5,7,-2)/1000),nrow=5,byrow=F)
sdE<-c(0.300000, 1.019467, 3, 2.030699, 2.328935)
X<-rep(1,N)
for (i in 1:5){
  confi<-rnorm(n=N,mean=Eobs%*%t(betaE[i,,drop=F]),sd=sdE[i])
  X<-cbind(X,confi)
}

## simulate 5 additional predictors of Y that aren't confounders and aren't used in the matching ##
predonly<-cbind(rnorm(N),rexp(N),runif(N),rnorm(N,sd=2.5))
betapred<-c(.04,-.02,.05,-.03)


## beta values and predictor forms ##
if (simtype==1){
  betaconf<-c(0.051400702, -0.0024289523,  0.0019607745, -0.0007199380,  0.0073785807,  0.043575231,-0.00013112471, -0.14723793)
  beta<-matrix(c(-1,betaconf,.008,.0001,.0027,.7/1000,.2/1000,.3/1000,betapred),ncol=1)  
  
  
  pform.cf<-cbind(X,X[,2]*X[,3],X[,4]^2,exp(X[,5])/(1+exp(X[,5])),Ecf,Ecf[,1]^2,Ecf[,1]*Ecf[,2],Ecf[,1]*Ecf[,2]*abs(X[,6]),Ecf[,1]*Ecf[,2]*abs(X[,5]),predonly)
  pform.obs<-cbind(X,X[,2]*X[,3],X[,4]^2,exp(X[,5])/(1+exp(X[,5])),Eobs,Eobs[,1]^2,Eobs[,1]*Eobs[,2],Eobs[,1]*Eobs[,2]*abs(X[,6]),Eobs[,1]*Eobs[,2]*abs(X[,5]),predonly)
  
} else if (simtype==2){ 
  betaconf<-c(0.051400702, -0.0024289523,  0.0019607745, -0.0007199380,  0.0073785807,  0.013575231,-0.013112471, -0.24723793)
  beta<-matrix(c(2,betaconf,.004,.0001,.007,.1/1000,betapred),ncol=1)  
  
  pform.cf<-cbind(X,X[,2]*X[,3],X[,4]^2,exp(X[,5])/(1+exp(X[,5])),Ecf,Ecf[,1]^2,Ecf[,1]*Ecf[,2],predonly)
  pform.obs<-cbind(X,X[,2]*X[,3],X[,4]^2,exp(X[,5])/(1+exp(X[,5])),Eobs,Eobs[,1]^2,Eobs[,1]*Eobs[,2],predonly)
  
  
} else if (simtype==3){
  betaconf<-c(0.001400702, -0.024289523,  0.049607745, -0.037199380,  0.023785807)
  beta<-matrix(c(5,betaconf,.01,.005,betapred),ncol=1)
  
  pform.cf<-cbind(X,Ecf,predonly)
  pform.obs<-cbind(X,Eobs,predonly)
  
}


## true causal effect for each unit ##
ice<-exp(pform.cf%*%beta)-exp(pform.obs%*%beta)

## true TEA for all units (tau) ##
trueTEAall<-sum(ice)

## omega values, tolerances for exact matching on E ##
ecut<-apply(Ecf,2,sd)*calp

xx<-matrix(ecut,nrow=nrow(Ecf),ncol=ncol(Ecf),byrow=T)

## caliper for mahalanobis distance matching on confounders ##
quan<-NULL
S_inv<-solve(cov(X[,-1]))
for (j in 1:N){
  temp<-matrix(rep(X[j,-1],N-1),nrow=N-1,byrow=T)
  quan<-c(quan,quantile(sqrt(rowSums((temp-X[-j,-1])%*%S_inv*(temp-X[-j,-1]))),.1))
}

## carry out matching ##
## the matching function outputs a dataframe with columns
## trtind: the indices of the units for which we found matches for their counterfactual pollutants and their confounders
## ctlind: the indices of the matched units for the unit given in the trtind column in the same row
## note that units for which multiple matches are found will be repeated in the trtind column,
## with each row having a different ctlind value corresponding to a different matched unit
orig.matchout.1<-match_rn_app_v2(trtobs=Eobs,trtcf=Ecf,confounders=X[,-1],trtdiff=xx,mdqtl=mean(quan))
orig.matchind.1<-as.data.frame(orig.matchout.1)
names(orig.matchind.1)<-c('trtind','ctlind')
## keep.1 is a vector of the indices of the retained units after trimming
keep.1<-unique(orig.matchind.1$trtind)

########################
## 2. RUN SIMULATIONS ##
########################

## initiate variables to save results for each method ##
match.1<-NULL
match.2<-NULL
match.3<-NULL

bart.1<-NULL

lm.1<-NULL
lm.2<-NULL

trueTEAretain<-NULL

for (g in 1:reps){
  ## for each set of simulated data, you want a different seed here in order to obtain different Y values ##
  set.seed(g-1)
  
  ## simulate observed counts of health events, associated with pollutants and confounders ##
  Yobs<-rpois(n=N,lambda=exp(pform.obs%*%beta))
  
  ######################################
  ## compute matching point estimates ##
  ######################################
  
  ## match.1=no bias correction ##
  ## impute counterfactuals ##
  EYcf<-tapply(Yobs[orig.matchind.1$ctlind],orig.matchind.1$trtind,mean)
  ## difference in imputed counterfactuals and observed number of deaths in each area ## 
  match.est.1<-sum(EYcf-Yobs[as.numeric(names(EYcf))])
  
  ## match.2=GAM bias correction ##
  regdat<-data.frame(Yobs,Eobs,X[,-1])
  names(regdat)<-c('y','pm25','o3',paste0('conf',1:5))
  fit.1<-bam(y~s(pm25)+s(o3)+s(conf1)+s(conf2)+s(conf3)+s(conf4)+s(conf5), family="poisson", data=regdat,discrete=T)
  ## for each CF level in matched data, use regression function to predict outcome with true confounders and with confounders of matches ##
  predXi<-data.frame(Ecf[orig.matchind.1$trtind,],X[orig.matchind.1$trtind,-1])
  predXj<-data.frame(Ecf[orig.matchind.1$trtind,],X[orig.matchind.1$ctlind,-1])
  names(predXi)<-c('pm25','o3',paste0('conf',1:5))
  names(predXj)<-c('pm25','o3',paste0('conf',1:5))
  ## take difference to be used for bias correction ##
  diffmeans<-predict(fit.1, predXi, type="response")-predict(fit.1, predXj, type="response")
  ## impute counterfactuals ##
  EYcf<-tapply(Yobs[orig.matchind.1$ctlind]+diffmeans,orig.matchind.1$trtind,mean)
  ## difference in imputed counterfactuals and observed number of deaths in each area ## 
  match.est.2<-sum(EYcf-Yobs[as.numeric(names(EYcf))])
  
  ## match.3=correct/oracle bias correction model ##
  if (simtype==1){
    regdat<-data.frame(Yobs,Eobs,X[,-1],exp(X[,5])/(1+exp(X[,5])),abs(X[,5:6]))
    names(regdat)<-c('y','pm25','o3',paste0('conf',1:8))
    fit.1<-glm(y~pm25+o3+I(pm25^2)+pm25*o3+conf1+conf2+conf3+conf4+conf5+
                 pm25*o3*conf7+pm25*o3*conf8+conf1*conf2+I(conf3^2)+conf6, family="poisson", data=regdat)
    ## for each CF level in matched data, use regression function to predict outcome with true confounders and with confounders of matches ##
    predXi<-data.frame(Ecf[orig.matchind.1$trtind,],X[orig.matchind.1$trtind,-1],exp(X[orig.matchind.1$trtind,5])/(1+exp(X[orig.matchind.1$trtind,5])),abs(X[orig.matchind.1$trtind,5:6]))
    predXj<-data.frame(Ecf[orig.matchind.1$trtind,],X[orig.matchind.1$ctlind,-1],exp(X[orig.matchind.1$ctlind,5])/(1+exp(X[orig.matchind.1$ctlind,5])),abs(X[orig.matchind.1$ctlind,5:6]))
    names(predXi)<-c('pm25','o3',paste0('conf',1:8))
    names(predXj)<-c('pm25','o3',paste0('conf',1:8))
    ## take difference to be used for bias correction ##
    diffmeans<-predict(fit.1, predXi, type="response")-predict(fit.1, predXj, type="response")
    ## impute counterfactuals ##
    EYcf<-tapply(Yobs[orig.matchind.1$ctlind]+diffmeans,orig.matchind.1$trtind,mean)
    ## difference in imputed counterfactuals and observed number of deaths in each area ## 
    match.est.3<-sum(EYcf-Yobs[as.numeric(names(EYcf))])
    
  } else if (simtype==2){
    ## match.3=correct model for simtype 2 bias correction ##
    regdat<-data.frame(Yobs,Eobs,X[,-1],exp(X[,5])/(1+exp(X[,5])))
    names(regdat)<-c('y','pm25','o3',paste0('conf',1:6))
    fit.1<-glm(y~pm25+o3+I(pm25^2)+pm25*o3+conf1+conf2+conf3+conf4+conf5+
                 conf1*conf2+I(conf3^2)+conf6, family="poisson", data=regdat)
    ## for each CF level in matched data, use regression function to predict outcome with true confounders and with confounders of matches ##
    predXi<-data.frame(Ecf[orig.matchind.1$trtind,],X[orig.matchind.1$trtind,-1],exp(X[orig.matchind.1$trtind,5])/(1+exp(X[orig.matchind.1$trtind,5])))
    predXj<-data.frame(Ecf[orig.matchind.1$trtind,],X[orig.matchind.1$ctlind,-1],exp(X[orig.matchind.1$ctlind,5])/(1+exp(X[orig.matchind.1$ctlind,5])))
    names(predXi)<-c('pm25','o3',paste0('conf',1:6))
    names(predXj)<-c('pm25','o3',paste0('conf',1:6))
    ## take difference to be used for bias correction ##
    diffmeans<-predict(fit.1, predXi, type="response")-predict(fit.1, predXj, type="response")
    ## impute counterfactuals ##
    EYcf<-tapply(Yobs[orig.matchind.1$ctlind]+diffmeans,orig.matchind.1$trtind,mean)
    ## difference in imputed counterfactuals and observed number of deaths in each area ## 
    match.est.3<-sum(EYcf-Yobs[as.numeric(names(EYcf))])
    
  } else if (simtype==3){
    ## match.3=linear model bias correction ##
    regdat<-data.frame(Yobs,Eobs,X[,-1])
    names(regdat)<-c('y','pm25','o3',paste0('conf',1:5))
    fit.1<-glm(y~pm25+o3+conf1+conf2+conf3+conf4+conf5, family="poisson", data=regdat)
    ## for each CF level in matched data, use regression function to predict outcome with true confounders and with confounders of matches ##
    predXi<-data.frame(Ecf[orig.matchind.1$trtind,],X[orig.matchind.1$trtind,-1])
    predXj<-data.frame(Ecf[orig.matchind.1$trtind,],X[orig.matchind.1$ctlind,-1])
    names(predXi)<-c('pm25','o3',paste0('conf',1:5))
    names(predXj)<-c('pm25','o3',paste0('conf',1:5))
    ## take difference to be used for bias correction ##
    diffmeans<-predict(fit.1, predXi, type="response")-predict(fit.1, predXj, type="response")
    ## impute counterfactuals ##
    EYcf<-tapply(Yobs[orig.matchind.1$ctlind]+diffmeans,orig.matchind.1$trtind,mean)
    ## difference in imputed counterfactuals and observed number of deaths in each area ## 
    match.est.3<-sum(EYcf-Yobs[as.numeric(names(EYcf))])
    
  }

  ######################
  ## do bootstrapping ##
  ######################
    
  match.1.bootest<-NULL
  match.2.bootest<-NULL
  match.3.bootest<-NULL
  
  match.boottru<-NULL

  bootind1<-list()
  bootind2<-list()
  for (j in 1:nboot){
    
    ## bootstrap from the trimmed sample ##
    bootind1<-c(bootind1,list(sample(x=keep.1,size=length(keep.1),replace=T)))
    
    ## bootstrap from the untrimmed sample ##
    bootind2<-c(bootind2,list(sample(x=1:N,size=N,replace=T)))
    
  }
    
  ## run the matching on the bootstrapped data as described in Section 3 of the paper ##
  bootmatch<-mapply(bootind1,bootind2,FUN=match_rn_boot6,MoreArgs=list(Eobs=Eobs,Ecf=Ecf,X=X[,-1],matchind=orig.matchind.1),SIMPLIFY=T)

  ##########################################################
  ## compute matching estimates for each bootstrap sample ##
  ##########################################################
    
  for (j in 1:nboot){
    bootmatchj<-bootmatch[[j]]
    
    ## match.1 ##
    EYcf<-tapply(Yobs[bootmatchj$ctlind],bootmatchj$unqind,mean)
    
    ## difference in imputed counterfactuals and observed number of deaths in each area ## 
    match.1.bootest<-c(match.1.bootest,sum(EYcf-Yobs[bootind1[[j]]]))
    
    ## match.2 ##
    regdat<-data.frame(Yobs[bootind2[[j]]],Eobs[bootind2[[j]],],X[bootind2[[j]],-1])
    names(regdat)<-c('y','pm25','o3',paste0('conf',1:5))
    fit.1<-bam(y~s(pm25)+s(o3)+s(conf1)+s(conf2)+s(conf3)+s(conf4)+s(conf5), family="poisson", data=regdat,discrete=T)
    ## for each CF level in matched data, use regression function to predict outcome with true confounders and with confounders of matches ##
    predXi<-data.frame(Ecf[bootmatchj$trtind,],X[bootmatchj$trtind,-1])
    predXj<-data.frame(Ecf[bootmatchj$trtind,],X[bootmatchj$ctlind,-1])
    names(predXi)<-c('pm25','o3',paste0('conf',1:5))
    names(predXj)<-c('pm25','o3',paste0('conf',1:5))
    ## take difference to be used for bias correction ##
    diffmeans<-predict(fit.1, predXi, type="response")-predict(fit.1, predXj, type="response")
    ## impute counterfactuals ##
    EYcf<-tapply(Yobs[bootmatchj$ctlind]+diffmeans,bootmatchj$unqind,mean)
    ## difference in imputed counterfactuals and observed number of deaths in each area ## 
    match.2.bootest<-c(match.2.bootest,sum(EYcf-Yobs[bootind1[[j]]]))
    
    if (simtype==1){
      ## match.3 ##
      regdat<-data.frame(Yobs[bootind2[[j]]],Eobs[bootind2[[j]],],X[bootind2[[j]],-1],exp(X[bootind2[[j]],5])/(1+exp(X[bootind2[[j]],5])),abs(X[bootind2[[j]],5:6]))
      names(regdat)<-c('y','pm25','o3',paste0('conf',1:8))
      fit.1<-glm(y~pm25+o3+I(pm25^2)+pm25*o3+conf1+conf2+conf3+conf4+conf5+
                   pm25*o3*conf7+pm25*o3*conf8+conf1*conf2+I(conf3^2)+conf6, family="poisson", data=regdat)
      ## for each CF level in matched data, use regression function to predict outcome with true confounders and with confounders of matches ##
      predXi<-data.frame(Ecf[bootmatchj$trtind,],X[bootmatchj$trtind,-1],exp(X[bootmatchj$trtind,5])/(1+exp(X[bootmatchj$trtind,5])),abs(X[bootmatchj$trtind,5:6]))
      predXj<-data.frame(Ecf[bootmatchj$trtind,],X[bootmatchj$ctlind,-1],exp(X[bootmatchj$ctlind,5])/(1+exp(X[bootmatchj$ctlind,5])),abs(X[bootmatchj$ctlind,5:6]))
      names(predXi)<-c('pm25','o3',paste0('conf',1:8))
      names(predXj)<-c('pm25','o3',paste0('conf',1:8))
      ## take difference to be used for bias correction ##
      diffmeans<-predict(fit.1, predXi, type="response")-predict(fit.1, predXj, type="response")
      ## impute counterfactuals ##
      EYcf<-tapply(Yobs[bootmatchj$ctlind]+diffmeans,bootmatchj$unqind,mean)
      ## difference in imputed counterfactuals and observed number of deaths in each area ## 
      match.3.bootest<-c(match.3.bootest,sum(EYcf-Yobs[bootind1[[j]]]))
      
    } else if (simtype==2){
      ## match.3 ##
      regdat<-data.frame(Yobs[bootind2[[j]]],Eobs[bootind2[[j]],],X[bootind2[[j]],-1],exp(X[bootind2[[j]],5])/(1+exp(X[bootind2[[j]],5])))
      names(regdat)<-c('y','pm25','o3',paste0('conf',1:6))
      fit.1<-glm(y~pm25+o3+I(pm25^2)+pm25*o3+conf1+conf2+conf3+conf4+conf5+
                   conf1*conf2+I(conf3^2)+conf6, family="poisson", data=regdat)
      ## for each CF level in matched data, use regression function to predict outcome with true confounders and with confounders of matches ##
      predXi<-data.frame(Ecf[bootmatchj$trtind,],X[bootmatchj$trtind,-1],exp(X[bootmatchj$trtind,5])/(1+exp(X[bootmatchj$trtind,5])))
      predXj<-data.frame(Ecf[bootmatchj$trtind,],X[bootmatchj$ctlind,-1],exp(X[bootmatchj$ctlind,5])/(1+exp(X[bootmatchj$ctlind,5])))
      names(predXi)<-c('pm25','o3',paste0('conf',1:6))
      names(predXj)<-c('pm25','o3',paste0('conf',1:6))
      ## take difference to be used for bias correction ##
      diffmeans<-predict(fit.1, predXi, type="response")-predict(fit.1, predXj, type="response")
      ## impute counterfactuals ##
      EYcf<-tapply(Yobs[bootmatchj$ctlind]+diffmeans,bootmatchj$unqind,mean)
      ## difference in imputed counterfactuals and observed number of deaths in each area ## 
      match.3.bootest<-c(match.3.bootest,sum(EYcf-Yobs[bootind1[[j]]]))
      
    } else if (simtype==3){
      regdat<-data.frame(Yobs[bootind2[[j]]],Eobs[bootind2[[j]],],X[bootind2[[j]],-1])
      names(regdat)<-c('y','pm25','o3',paste0('conf',1:5))
      fit.1<-glm(y~pm25+o3+conf1+conf2+conf3+conf4+conf5, family="poisson", data=regdat)
      ## for each CF level in matched data, use regression function to predict outcome with true confounders and with confounders of matches ##
      predXi<-data.frame(Ecf[bootmatchj$trtind,],X[bootmatchj$trtind,-1])
      predXj<-data.frame(Ecf[bootmatchj$trtind,],X[bootmatchj$ctlind,-1])
      names(predXi)<-c('pm25','o3',paste0('conf',1:5))
      names(predXj)<-c('pm25','o3',paste0('conf',1:5))
      ## take difference to be used for bias correction ##
      diffmeans<-predict(fit.1, predXi, type="response")-predict(fit.1, predXj, type="response")
      ## impute counterfactuals ##
      EYcf<-tapply(Yobs[bootmatchj$ctlind]+diffmeans,bootmatchj$unqind,mean)
      ## difference in imputed counterfactuals and observed number of deaths in each area ## 
      match.3.bootest<-c(match.3.bootest,sum(EYcf-Yobs[bootind1[[j]]]))
      
    }
    
    ## save the truth for each bootstrap sample ##
    match.boottru<-c(match.boottru,sum(ice[bootind1[[j]]]))
  }
  bootmatch<-NULL
  bootmatchj<-NULL

  ## save the matching results and print to console ##
  match.1<-rbind(match.1,c(match.est.1,quantile(match.1.bootest,.025),quantile(match.1.bootest,.975)))
  match.2<-rbind(match.2,c(match.est.2,quantile(match.2.bootest,.025),quantile(match.2.bootest,.975)))
  match.3<-rbind(match.3,c(match.est.3,quantile(match.3.bootest,.025),quantile(match.3.bootest,.975)))

  print(paste0('True TEA: ',sum(ice[keep.1])))
  print(paste0('Matching with no bias correction: ',match.est.1,' (',quantile(match.1.bootest,.025),',',quantile(match.1.bootest,.975),')'))
  print(paste0('Matching with GAM bias correction: ',match.est.2,' (',quantile(match.2.bootest,.025),',',quantile(match.2.bootest,.975),')'))
  print(paste0('Matching with oracle bias correction: ',match.est.3,' (',quantile(match.3.bootest,.025),',',quantile(match.3.bootest,.975),')'))
  
  
  ##########
  ## BART ##
  ##########
  
  ## fit the BART model to the observed data ##
  ## use BART to predict counterfactuals ##
  bartfit<-bart(x.train=cbind(Eobs,X[,-1]),y.train=as.numeric(Yobs),x.test=cbind(Ecf,X[,-1]),verbose = F)
  
  ## use BART to estimate the TEA (note that we are only using the units retained after trimming for estimation) ##
  bartest_tr<-sum(bartfit$yhat.test.mean[keep.1]-Yobs[keep.1])
  
  ## compute the credible interval ##
  tau_pd_tr<-NULL
  for (j in 1:1000){
    tau_pd_tr<-c(tau_pd_tr,sum(rnorm(n=length(keep.1),mean=bartfit$yhat.test[j,keep.1],sd=bartfit$sigma[j])-Yobs[keep.1]))
  }
  
  ## save the BART results and print to console ##
  bart.1<-rbind(bart.1, c(bartest_tr,quantile(tau_pd_tr,c(.025,.975))))
  print(paste0('BART: ',bartest_tr,' (',quantile(tau_pd_tr,.025),',',quantile(tau_pd_tr,.975),')'))
  
  ########################
  ## Poisson Regression ##
  ########################
   
  ## lm.1=linear terms only ##
  regdat<-data.frame(Yobs,Eobs,X[,-1])
  names(regdat)<-c('y','pm25','o3',paste0('conf',1:5))
  fit.3<-glm(y~pm25+o3+conf1+conf2+conf3+conf4+conf5, family="poisson", data=regdat)
  ## predict at counterfactual pollution levels for all units ##
  predXi<-data.frame(Ecf[keep.1,],X[keep.1,-1])
  names(predXi)<-c('pm25','o3',paste0('conf',1:5))
  EYcf<-predict(fit.3, predXi, type="response")
  ice_hat.3<-EYcf-Yobs[keep.1]
  ## compute var-cov matrix for the predictions so that we can get the variance of the their sum ##
  vdat<-data.frame(Yobs[keep.1],Ecf[keep.1,],X[keep.1,-1])
  names(vdat)<-c('y','pm25','o3',paste0('conf',1:5))
  foo<-model.matrix(y~pm25+o3+conf1+conf2+conf3+conf4+conf5,data=vdat)
  pred_vcov<-(foo*matrix(rep(EYcf,ncol(foo)),nrow=length(keep.1),ncol=ncol(foo),byrow=F))%*%vcov(fit.3)%*%
            t(foo*matrix(rep(EYcf,ncol(foo)),nrow=length(keep.1),ncol=ncol(foo),byrow=F))
  
  ## save results and print to console ##
  lm.1<-rbind(lm.1,c(sum(ice_hat.3),
                           sum(ice_hat.3)-1.96*sqrt(matrix(1,ncol=length(keep.1),nrow=1)%*%pred_vcov%*%matrix(1,ncol=1,nrow=length(keep.1))),
                           sum(ice_hat.3)+1.96*sqrt(matrix(1,ncol=length(keep.1),nrow=1)%*%pred_vcov%*%matrix(1,ncol=1,nrow=length(keep.1)))))
  print(paste0('PR with linear terms: ',lm.1[nrow(lm.1),1],' (',lm.1[nrow(lm.1),2],',',lm.1[nrow(lm.1),3],')'))
  
  ## lm.2=correct model for simtypes 1 and 2 ##
  if (simtype==1){
    regdat<-data.frame(Yobs,Eobs,X[,-1],exp(X[,5])/(1+exp(X[,5])),abs(X[,5:6]))
    names(regdat)<-c('y','pm25','o3',paste0('conf',1:8))
    fit.3<-glm(y~pm25+o3+I(pm25^2)+pm25*o3+conf1+conf2+conf3+conf4+conf5+
                 pm25*o3*conf7+pm25*o3*conf8+conf1*conf2+I(conf3^2)+conf6, family="poisson", data=regdat)
    ## predict at counterfactual pollution levels for all units ##
    predXi<-data.frame(Ecf[keep.1,],X[keep.1,-1],exp(X[keep.1,5])/(1+exp(X[keep.1,5])),abs(X[keep.1,5:6]))
    names(predXi)<-c('pm25','o3',paste0('conf',1:8))
    EYcf<-predict(fit.3, predXi, type="response")
    ice_hat.3<-EYcf-Yobs[keep.1]
    ## compute var-cov matrix for the predictions so that we can get the variance of the their sum ##
    vdat<-data.frame(Yobs[keep.1],Ecf[keep.1,],X[keep.1,-1],exp(X[keep.1,5])/(1+exp(X[keep.1,5])),abs(X[keep.1,5:6]))
    names(vdat)<-c('y','pm25','o3',paste0('conf',1:8))
    foo<-model.matrix(y~pm25+o3+I(pm25^2)+pm25*o3+conf1+conf2+conf3+conf4+conf5+
                        pm25*o3*conf7+pm25*o3*conf8+conf1*conf2+I(conf3^2)+conf6,data=vdat)
    pred_vcov<-(foo*matrix(rep(EYcf,ncol(foo)),nrow=length(keep.1),ncol=ncol(foo),byrow=F))%*%vcov(fit.3)%*%
      t(foo*matrix(rep(EYcf,ncol(foo)),nrow=length(keep.1),ncol=ncol(foo),byrow=F))
    
    ## save results and print to console ##
    lm.2<-rbind(lm.2,c(sum(ice_hat.3),
                       sum(ice_hat.3)-1.96*sqrt(matrix(1,ncol=length(keep.1),nrow=1)%*%pred_vcov%*%matrix(1,ncol=1,nrow=length(keep.1))),
                       sum(ice_hat.3)+1.96*sqrt(matrix(1,ncol=length(keep.1),nrow=1)%*%pred_vcov%*%matrix(1,ncol=1,nrow=length(keep.1)))))
    print(paste0('PR correctly specified: ',lm.2[nrow(lm.2),1],' (',lm.2[nrow(lm.2),2],',',lm.2[nrow(lm.2),3],')'))
    
  } else if (simtype==2){
    ## lm.2=correct model for simtype 2 ##
    regdat<-data.frame(Yobs,Eobs,X[,-1],exp(X[,5])/(1+exp(X[,5])))
    names(regdat)<-c('y','pm25','o3',paste0('conf',1:6))
    fit.3<-glm(y~pm25+o3+I(pm25^2)+pm25*o3+conf1+conf2+conf3+conf4+conf5+
                 conf1*conf2+I(conf3^2)+conf6, family="poisson", data=regdat)
    ## predict at counterfactual pollution levels for all units ##
    predXi<-data.frame(Ecf[keep.1,],X[keep.1,-1],exp(X[keep.1,5])/(1+exp(X[keep.1,5])))
    names(predXi)<-c('pm25','o3',paste0('conf',1:6))
    EYcf<-predict(fit.3, predXi, type="response")
    ice_hat.3<-EYcf-Yobs[keep.1]
    ## compute var-cov matrix for the predictions so that we can get the variance of the their sum ##
    vdat<-data.frame(Yobs[keep.1],Ecf[keep.1,],X[keep.1,-1],exp(X[keep.1,5])/(1+exp(X[keep.1,5])))
    names(vdat)<-c('y','pm25','o3',paste0('conf',1:6))
    foo<-model.matrix(y~pm25+o3+I(pm25^2)+pm25*o3+conf1+conf2+conf3+conf4+conf5+
                        conf1*conf2+I(conf3^2)+conf6,data=vdat)
    pred_vcov<-(foo*matrix(rep(EYcf,ncol(foo)),nrow=length(keep.1),ncol=ncol(foo),byrow=F))%*%vcov(fit.3)%*%
      t(foo*matrix(rep(EYcf,ncol(foo)),nrow=length(keep.1),ncol=ncol(foo),byrow=F))
    
    ## save results and print to console ##
    lm.2<-rbind(lm.2,c(sum(ice_hat.3),
                       sum(ice_hat.3)-1.96*sqrt(matrix(1,ncol=length(keep.1),nrow=1)%*%pred_vcov%*%matrix(1,ncol=1,nrow=length(keep.1))),
                       sum(ice_hat.3)+1.96*sqrt(matrix(1,ncol=length(keep.1),nrow=1)%*%pred_vcov%*%matrix(1,ncol=1,nrow=length(keep.1)))))
    print(paste0('PR correctly specified: ',lm.2[nrow(lm.2),1],' (',lm.2[nrow(lm.2),2],',',lm.2[nrow(lm.2),3],')'))
    
  }

  ## true TEA for only units retained after trimming ##
  trueTEAretain<-c(trueTEAretain,sum(ice[keep.1]))

}

## if desired, use the code below to save the results to an external .RData file ##

match.1<-data.frame('Match 1',match.1)
match.2<-data.frame('Match 2',match.2)
match.3<-data.frame('Match 3',match.3)

bart.1<-data.frame('BART',bart.1)

lm.1<-data.frame('PR 1',lm.1)

names(match.1)<-names(match.2)<-names(match.3)<-c('method','est','ll','ul')

names(bart.1)<-c('method','est','ll','ul')

names(lm.1)<-c('method','est','ll','ul')
if (length(lm.2)>0){
  lm.2<-data.frame('PR 2',lm.2)
  names(lm.2)<-c('method','est','ll','ul')
}

save(trueTEAall,trueTEAretain,match.1,match.2,match.3,bart.1,lm.1,lm.2,
     file=paste0('sim',simtype,'_tol',calp*100,'_results.RData'))

