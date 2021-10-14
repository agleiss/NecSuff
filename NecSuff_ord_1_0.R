##################################################################
#
# NecSuff_ord
# ===========
#
# R-function for computing degrees of necessity and sufficiency
# for an ordinal outcome
#
# Gleiss, A., Henderson, R., Schemper, M., Degrees of necessity 
# and of sufficiency: further results and extensions, with an 
# application to covid-19 mortality in Austria
# accepted by Statistics in Medicine
#
# Author:  Andreas Gleiss
# Version: 1.0
# Date:    14 Oct 2020
#
# Arguments:
# ==========
# 
# pred 			name of variable containing predictions
#           (no. of columns = no. of outcome categories,
#            assumes that ordinal outcome is sorted from "best"
#            to "worst" outcome category)
#
##################################################################

NecSuff_ord<-function(pred){
  nc<-ncol(pred)
  nsi<-matrix(NA,nc,6)
  wgt<-apply(pred,2,sum)[2:nc]/sum(pred[,2:nc])
  nsi[,6]<-c(wgt,0)
  colnames(nsi)<-c("DN1","DS1","DN2","DS2","EV","weight")
  rownames(nsi)<-c(as.character(2:nc),"w.sum")
  for (i in 2:nc) {
    if(i<nc)
      predi<-apply(pred[,i:nc],1,sum)
    else
      predi<-pred[,nc]
    nsi[(i-1),1:5]<-NecSuff(predi, print=0)
  }
  nsi[nc,1:5]<-wgt%*%nsi[1:(nc-1),1:5]
  
  print(round(nsi,3))
}

xx<-rnorm(1000)
y<-5+0.5*xx+0.5*rnorm(1000)
yy<-(y>4)+(y>5)+(y>6) # vgl. Agresti
plot(y~xx, col=yy+1)
yyf<-factor(yy, ordered=T)
model<-polr(yyf~xx, Hess = TRUE)
summary(model)

pred<-predict(model, xx, type = "probs")
NecSuff_ord(pred)

#goeg$ord_<-factor(goeg$ord, ordered=T)      # best to worst
#model<-polr(goeg$ord_~goeg$alter65, Hess = TRUE) 
#summary(model) 
#pred<-predict(model, goeg$alter65, type = "probs")
#NecSuff_ord(pred)


