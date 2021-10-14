##################################################################
#
# NecSuff
# =======
#
# R-function for computing degrees of necessity and sufficiency
# for a dichotomous outcome
#
# Gleiss, A. & Schemper, M. Quantifying degrees of necessity and 
# sufficiency in cause-effect relationships with dichotomous and 
# survival outcome,
# Statistics in Medicine 2019; 38:4733-4748.
#
# Author:  Andreas Gleiss
# Version: 1.2 (output for use in NecSuff_ord; print option)
# Date:    14 Oct 2020
#
# Arguments:
# ==========
# 
# pred 			name of variable containing predictions
# print     print output if =1 (default)
#
##################################################################

NecSuff<-function(pred, print=1){
  p_bar<-mean(pred)
  smaller<-(pred<p_bar)
  larger<-(pred>p_bar)
  DN1<-sqrt(mean(((p_bar-pred[smaller])/p_bar)^2))
  DS1<-sqrt(mean(((pred[larger]-p_bar)/(1-p_bar))^2))
  DN2<-mean((p_bar-pred[smaller])/p_bar)
  DS2<-mean((pred[larger]-p_bar)/(1-p_bar))
  
  w_dn<-p_bar/(1-p_bar) * sum(smaller)/length(pred)
  w_ds<-(1-p_bar)/p_bar * sum(larger)/length(pred)
  EV<-w_dn*DN1^2 + w_ds*DS1^2
  
  if(print) {
    cat('\nest.P(D) =',round(p_bar,3))
    cat('\nDN1 =',round(DN1,3),', DS1 =',round(DS1,3))
    cat('\nDN2 =',round(DN2,3),', DS1 =',round(DS2,3))
    cat('\nEV =',round(EV,3))
  }
  else
    return(c(DN1,DS1,DN2,DS2,EV))
  
  }

xx<-rnorm(1000)
pp<-1/(1+exp(-5-log(5)*xx))
yy<-1*(runif(1000)<=pp)
plot(yy~xx)
points(cbind(xx,pp),col="red")
res.glm<-glm(yy~xx, family=binomial)
summary(res.glm)
pred<-predict(res.glm, type="response")
NecSuff(pred)

# Swedish cohort (Nilsson et al, 2001)
xx<-c(rep(0,36),rep(1,177),rep(0,8120),rep(1,4331))
yy<-c(rep(1,36),rep(1,177),rep(0,8120),rep(0,4331))
table(yy,xx)
res.glm<-glm(yy~xx, family=binomial)
summary(res.glm)
pred<-predict(res.glm, type="response")
NecSuff(pred)
