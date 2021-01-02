## write a little function to (1) select treatment matches and (2) among treatment matches, select confounder matches ##
match_rn_boot6<-function(ind1,ind2,Eobs,Ecf,X,matchind){
  
  m.ind<-NULL
  for (oo in 1:length(ind1)){
    jj<-ind1[oo]
    foo<-matchind$ctlind[which(matchind$trtind==jj)]
    if (sum(foo %in% ind2)>0){
      m.ind<-rbind(m.ind,data.frame('unqind'=oo,'trtind'=jj,'ctlind'=ind2[which(ind2 %in% foo)]))
    }else{
      m.ind<-rbind(m.ind,data.frame('unqind'=oo,'trtind'=jj,'ctlind'=NA))
    }
  }
  
  ## for units with no matches, require them to have one match ##
  ## to find this match, get 20 closest units on treatment and choose the one that is closest in terms of confounders ##
  Ebootctl<-Eobs[ind2,]
  Xbootctl<-X[ind2,]
  
  S_inv_E<-solve(cov(Ebootctl))
  S_inv_X<-solve(cov(Xbootctl))
  
  for (i in which(is.na(m.ind$ctlind))){
    jj<-m.ind$trtind[i]
    
    temp<-matrix(as.numeric(Ecf[jj,]),nrow=nrow(Ebootctl),ncol=ncol(Ebootctl),byrow=T)
    MD_E<-sqrt(rowSums(as.matrix(temp-Ebootctl)%*%S_inv_E*(as.matrix(temp-Ebootctl))))
    ## 20 smallest MD's ##
    small_mde<-order(MD_E)[1:20]
    
    temp<-matrix(as.numeric(X[jj,]),nrow=length(small_mde),ncol=ncol(Xbootctl),byrow=T)
    MD_i<-sqrt(rowSums(as.matrix(temp-Xbootctl[small_mde,])%*%S_inv_X*(as.matrix(temp-Xbootctl[small_mde,]))))
    
    m.ind$ctlind[i]<-ind2[small_mde[which.min(MD_i)]]
  }

  return(list(m.ind))
  
}
