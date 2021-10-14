##################################################################
#
# NecSuff_nom0
# ============
#
# R-function for computing degrees of necessity and sufficiency
# for an nominal outcome with reference category
#
# Gleiss, A., Henderson, R., Schemper, M., Degrees of necessity 
# and of sufficiency: further results and extensions, with an 
# application to covid-19 mortality in Austria
# accepted by Statistics in Medicine
#
# Author:  Andreas Gleiss
# Version: 1.0
# Date:    15 Oct 2020
#
# Arguments:
# ==========
# 
# pred 			name of variable containing predictions
#           (no. of columns = no. of outcome categories)
# refcat    reference category (default = column 1) 
#
##################################################################

NecSuff_nom0<-function(pred, refcat=1){
  nc<-ncol(pred)
  nsi<-matrix(NA,nc,6)
  wgt<-apply(pred,2,sum)[2:nc]/sum(pred[,2:nc])
  nsi[,6]<-c(wgt,0)
  colnames(nsi)<-c("DN1","DS1","DN2","DS2","EV","weight")
  rownames(nsi)<-c(as.character(2:nc),"w.sum")
  for (i in 1:nc) {
    if(i!=refcat) {
      predi<-pred[,i]/apply(pred[,c(i,refcat)],1,sum)
      nsi[(i-1),1:5]<-NecSuff(predi, print=0)
    }
  }
  nsi[nc,1:5]<-wgt%*%nsi[1:(nc-1),1:5]
  
  print(round(nsi,3))
}

#library(nnet)

#goeg$nom0<-relevel(factor(goeg$ord), ref="0")
#model<-multinom(goeg$nom0~goeg$alter65, Hess = TRUE) 
#summary(model) 
#pred<-predict(model, goeg$alter65, type = "probs")
#NecSuff_nom0(pred)


