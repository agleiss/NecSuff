##################################################################
#
# NecSuff_CR
# ==========
#
# R-function for computing degrees of necessity and sufficiency
# for a survival outcome (variant 1) with competing risk based on
# Fine & Gray model, together with the direct V measure and the
# indirect V_W measure as proportion of explained variation
#
# Gleiss, A., Gnant, M., Schemper, M., Explained variation and 
# degrees of necessity and of sufficiency for competing risks 
# survival data. Submitted to SMMR.
#
# Author:  Andreas Gleiss
# Version: 1.0 
# Date:    17 Mar 2022
#
# Arguments:
# ==========
# 
# time      time variable
# status 	  status variable (censored values have status 0)
# x         independent variable   
# failcode  code for event type of interest (default=1)
#
##################################################################

NecSuff_CR <- function (time, status, x, failcode=1) 
{
  f.Mt.ind <- function(tempo, tutti.tempi,   
                       Cti, Ctj,  
                       Ctix, pp,
                       tempi.evento, status, failcode, num.sogg, 
                       Gtj, Gti) {
    
    Ctj1 <- unique(Ctj[tempi.evento == tempo])
    Ctjx<-pp[tempi.evento == tempo,2:(num.sogg+1)]
    
    qij <- Gtj[tempi.evento == tempo] / Gti
    qij[!(tutti.tempi <= tempo & status!=0 & status!=failcode)]<-1

    # for V_ind
    ris.V_num <- qij*Ctjx*(1-Ctjx)
    ris.V_den <- (Ctj1*(1-Ctj1))
    
    # for V_w_ind
    ris.V_W <-qij*Ctjx*(1-Ctjx) / (Ctj1*(1-Ctj1))
    
    # for DN
    rd <- qij*((Ctj1 - Ctjx) / Ctj1)^2
    ris.DN <- rd[(Ctjx<Ctj1) & (Ctj1>0)]

    # for DS
    rs <- qij*((Ctjx - Ctj1) / (1-Ctj1))^2
    ris.DS <- rs[(Ctjx>Ctj1) & (Ctj1<1)]
    
    return(c(ris.V_den, sum(ris.V_num)/sum(qij), 
             1-sum(ris.V_W)/sum(qij), 
             sqrt(sum(ris.DN)/sum((Ctjx<Ctj1)*qij)), 
             sqrt(sum(ris.DS)/sum((Ctjx>Ctj1)*qij)))) # sum gives i-loop in SAS
  }

  
  #function which calculates  M(t_(j)) for a given time t_(j)
  f.Mt <- function(tempo, tutti.tempi, 
                   Cti, tempi.evento, 
                   Ctj, status, failcode, num.sogg, 
                   Gtj, Gti)
  {
    Ctj1 <- unique(Ctj[tempi.evento == tempo])
    
    qij <- Gtj[tempi.evento == tempo] / Gti
    qij[!(tutti.tempi <= tempo & status!=0 & status!=failcode)]<-1

    primo <- rep(Ctj1, num.sogg)
    primo[tutti.tempi <= tempo] <- 0
    terzo <- (status==0) * (Ctj1*(1 - Ctj1)/(1-Cti) + (1-Ctj1)*(1 - (1-Ctj1)/(1-Cti)))
    terzo[tutti.tempi > tempo] <- 0
    terzo[is.na(terzo)] <- 0
    secondo <- (1-Ctj1) * (status==failcode)
    secondo[tutti.tempi > tempo] <- 0
    quatro <- Ctj1 * (status!=0 & status!=failcode) * qij
    quatro[tutti.tempi > tempo] <- 0
    ris <- primo + secondo + terzo + quatro
    return(sum(ris)/sum(qij)) # sum gives i-loop
  }
  
  #function which calculates M(t_(j)|x) for a given time t_(j)
  f.Mt.fg <- function(tempo, tutti.tempi, 
                      Ctix, tempi.evento, 
                      pp, status, failcode, num.sogg, 
                      Gtj, Gti)
  {
    Ctjx<-pp[tempi.evento == tempo,2:(num.sogg+1)]
    
    qij <- Gtj[tempi.evento == tempo] / Gti
    qij[!(tutti.tempi <= tempo & status!=0 & status!=failcode)]<-1
    
    primo <- Ctjx
    primo[tutti.tempi <= tempo] <- 0
    terzo <- (status==0) * (Ctjx*(1 - Ctjx)/(1-Ctix) + (1-Ctjx)*(1 - (1-Ctjx)/(1-Ctix)))
    terzo[tutti.tempi > tempo] <- 0
    terzo[is.na(terzo)] <- 0
    secondo <- (1-Ctjx) * (status==failcode)
    secondo[tutti.tempi > tempo] <- 0
    quatro <- Ctjx * (status!=0 & status!=failcode) * qij
    quatro[tutti.tempi > tempo] <- 0
    ris <- primo + secondo + terzo + quatro
    return(sum(ris)/sum(qij)) # sum gives i-loop
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
  dat<-na.omit(cbind(time,status,x)) # 1980 Zeilen
  time<-dat[,1]
  status<-dat[,2]
  x<-dat[,3]
  
  num.sogg <- length(time)
  cif.unc <- cuminc(time,status)
  tempi.eventi <- sort(unique(time[status==1])) # different EoI times
  pos.cif <- apply(as.matrix(time), 1, f.assegna.surv, tempi.eventi)
  unique.time.est<-timepoints(cif.unc,tempi.eventi)[[failcode]][1,] 

  surv.tot.cif.unc <- unique.time.est[pos.cif] 
  ind.censura <- as.numeric(!as.logical(status))
  cif.tj <- unique.time.est 
  numero.eventi <- as.vector(table(sort(time[status==failcode])))
  
  fgfit<-crr(time,status, x, failcode=failcode)
  pp <- predict(fgfit, x) # (1st column: tempi.eventi)
  cif.tot.fg <- pp[cbind(pos.cif,2:(num.sogg+1))]

  Gkm <- survfit(Surv(time, ind.censura) ~ 1) # reverse KM
  tempi.censure <- Gkm$time[Gkm$n.event != 0]
  if(!length(tempi.censure)) 
    cens.tot.km <- rep(1, length(tempi.eventi))
  else {
    pos.surv.censure <- apply(as.matrix(tempi.eventi), 1, f.assegna.surv, tempi.censure)
    cens.tot.km <- (Gkm$surv[Gkm$n.event != 0])[pos.surv.censure]
    cens.tot.km[tempi.eventi < min(Gkm$time[Gkm$n.event != 0])] <- 1
  } # G_hat(t(j))
  
  Mt <- apply(as.matrix(tempi.eventi), 1, f.Mt, time, 
              surv.tot.cif.unc, tempi.eventi, 
              cif.tj, status, failcode, num.sogg,
              cens.tot.km, cens.tot.km[pos.cif]) # apply gives j-loop
  Mtx <- apply(as.matrix(tempi.eventi), 1, f.Mt.fg, time, 
               cif.tot.fg, tempi.eventi, 
               pp, status, failcode, num.sogg, 
               cens.tot.km, cens.tot.km[pos.cif]) # apply gives j-loop

  pesi <- numero.eventi/cens.tot.km
  peso.tot <- sum(pesi)
  
  D <- sum(Mt * pesi)/peso.tot
  Dx <- sum(Mtx * pesi)/peso.tot
  V <- (D - Dx)/D
  if(any(Mt == 0)) {
    Vw <- sum((Mt[Mt != 0] - Mtx[Mt != 0])/Mt[Mt != 0] * pesi[Mt !=0])/sum(pesi[Mt != 0])}
  else Vw <- sum((Mt - Mtx)/Mt * pesi)/peso.tot

  
  # aus NecSuff_surv für indir. V_W, DN und DS
  Mt.ind <- apply(as.matrix(tempi.eventi), 1, f.Mt.ind, time,
                  surv.tot.cif.unc, cif.tj,  
                  cif.tot.fg, pp,
                  tempi.eventi, status, failcode, num.sogg, 
                  cens.tot.km, cens.tot.km[pos.cif]) # entspricht j-Schleife in SAS
  D_ind <- sum(Mt.ind[1,] * pesi)/peso.tot
  Dx_ind <- sum(Mt.ind[2,] * pesi)/peso.tot
  V_ind <- (D_ind - Dx_ind)/D_ind
  Vw_ind <- sum(Mt.ind[3,] * pesi)/peso.tot
  DN1 <- sum(Mt.ind[4,] * pesi)/peso.tot
  DS1 <- sum(Mt.ind[5,] * pesi)/peso.tot

  return(list(V = V, #V_ind = V_ind, Vw = Vw, 
              Vw_ind = Vw_ind, DN = DN1, DS = DS1))
}


  

### example

library(cmprsk)

#cif.unc<-cuminc(abcsg$tt_dist, abcsg$stat_dist)
#  plot(cif.unc)
#fg.tusize<-crr(abcsg$tt_dist, abcsg$stat_dist, abcsg$log_TumorSizeN)
#  summary(fg.tusize)
#NecSuff_CR(abcsg$tt_dist, abcsg$stat_dist, abcsg$log_TumorSizeN)

