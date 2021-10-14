##################################################################
#
# NecSuff_Surv
# ============
#
# R-function for computing degrees of necessity and sufficiency
# for a survival outcome (variant 1) together with the indirect
# V_W measure as proportion of explained variation
#
# Gleiss, A. & Schemper, M. Quantifying degrees of necessity and 
# sufficiency in cause-effect relationships with dichotomous and 
# survival outcome,
# Statistics in Medicine 2019; 38:4733-4748.
#
# Schemper M. Predictive accuracy and explained variation. 
# Statist Med. 2003;22(14):2299-2308.
#
# Author:  Andreas Gleiss
# Version: 1.0 
# Date:    17 Feb 2021
#
# Arguments:
# ==========
# 
# coxfit 	  result of cph()
#
##################################################################

NecSuff_Surv <- function (coxfit) 
{
  f.Mt.both <- function(tempo, tutti.tempi,   
                        Stj,
                        Stj0, lin.pred, 
                        tempi.evento,ind.censura, 
                        num.sogg) {
    KM <- unique(Stj[tempi.evento == tempo])
    
    Stj00 <- unique(Stj0[tempi.evento == tempo])
    Cox <- Stj00^exp(lin.pred)
    
    # for V_w_indir
    ris <- 1 - Cox*(1-Cox) / (KM*(1-KM))
    
    # for DN
    ris.DN <- (((Cox[(Cox>KM) & (KM<1)]-KM)/(1-KM))^2)

    # for DS
    ris.DS <- (((KM-Cox[(Cox<KM) & (KM>0)])/KM)^2)
    
    return(c(sum(ris)/num.sogg, 
             sqrt(sum(ris.DN)/sum(Cox>KM)), 
             sqrt(sum(ris.DS)/sum(Cox<KM)))) # sum entspricht i-Schleife in SAS
  }
  
  f.assegna.surv <- function(tempo, tempi.eventi) {
    if (any(tempo == tempi.eventi)) {
      pos <- (c(1:length(tempi.eventi)) * as.numeric(tempo == 
                                                       tempi.eventi))
      pos <- pos[pos != 0]
    }
    else {
      tmp <- (tempo - tempi.eventi)
      if (all(tmp < 0)) 
        pos <- NA
      else {
        tmp[tmp < 0] <- Inf
        pos <- order(tmp)[1]
      }
    }
    return(pos)
  }
  
  # begin of function
  tsurv <- as.numeric(coxfit$y[, 1])
  surv <- as.numeric(coxfit$y[, 2])
  num.sogg <- length(tsurv)
  km <- survfit(Surv(tsurv, surv) ~ 1)
  tempi.eventi <- km$time[km$n.event != 0]
  pos.surv <- apply(as.matrix(tsurv), 1, f.assegna.surv, tempi.eventi)
  surv.tot.km <- (km$surv[km$n.event != 0])[pos.surv]
  surv.tot.km[is.na(surv.tot.km)] <- 1
  ind.censura <- as.numeric(!as.logical(surv))
  surv.tj <- km$surv[km$n.event != 0]
  
  numero.eventi <- km$n.event[km$n.event != 0]
  surv0.tj.cox <- coxfit$surv[which(coxfit$time%in%tempi.eventi)]
  surv0.tot.cox <- (coxfit$surv[which(coxfit$time%in%tempi.eventi)])[pos.surv]
  surv.tot.cox <- surv0.tot.cox^exp(coxfit$linear.predictors)
  
  Mt.both <- apply(as.matrix(tempi.eventi), 1, f.Mt.both, 
                   surv.tot.km, surv.tj,
                   surv0.tj.cox, coxfit$linear.predictors, 
                   tempi.eventi, ind.censura, 
                   num.sogg) # entspricht j-Schleife in SAS

  Gkm <- survfit(Surv(tsurv, ind.censura) ~ 1)
  tempi.censure <- Gkm$time[Gkm$n.event != 0]
  if (!length(tempi.censure)) 
    cens.tot.km <- rep(1, length(tempi.eventi))
  else {
    pos.surv.censure <- apply(as.matrix(tempi.eventi), 1, 
                              f.assegna.surv, tempi.censure)
    cens.tot.km <- (Gkm$surv[Gkm$n.event != 0])[pos.surv.censure]
    cens.tot.km[tempi.eventi < min(Gkm$time[Gkm$n.event != 
                                              0])] <- 1
  }
  pesi <- numero.eventi/cens.tot.km
  peso.tot <- sum(pesi)
  
  Vw_indir <- sum(Mt.both[1,] * pesi)/peso.tot
  
  DN1 <- sum(Mt.both[2,] * pesi)/peso.tot
  
  DS1 <- sum(Mt.both[3,] * pesi)/peso.tot
  
  return(list(Model = coxfit$call, Vw = Vw_indir, DN = DN1, DS = DS1))
}

#library(rms)
#surv<-Surv(pbcneu$zeit,pbcneu$status)
#plot(survfit(surv~pbcneu$edema))
#coxfit<-cph(surv~age+edema, y=TRUE, surv=TRUE, method="breslow", type="kaplan-meier",
#            data=pbcneu)
#ens<-NecSuff_Surv(coxfit)
#print(ens)
