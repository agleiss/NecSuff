##################################################################
#
# NecSuff_CR
# ==========
#
# R-function for computing degrees of necessity and sufficiency
# for a survival outcome (variant 1) with competing risk based on
# Fine & Gray model, together with the direct and indirect V 
# measure and direct and indirect V_W measure as proportion of 
# explained variation
#
# Gleiss, A., Gnant, M., Schemper, M., Explained variation and 
# degrees of necessity and of sufficiency for competing risks 
# survival data. Submitted to Biometrical Journal.
#
# Author:  Andreas Gleiss
# Version: 1.1 (removed q and corrected p_ij)
# Date:    03 Oct 2023
#
# Arguments:
# ==========
# 
# time      time variable
# status 	  status variable (censored values have status 0)
# x         independent variable   
# failcode  code for event type of interest (default=1)
# CEcode    code for competing event (default=2)
#
##################################################################

NecSuff_CR <- function (time, status, x, failcode=1, CEcode=2) 
{
  f.Mt.ind <- function(tempo, tutti.tempi,   
                       Cti, Ctj,  
                       Ctix, pp,
                       tempi.evento, status, failcode, num.sogg, 
                       Gtj, Gti) {
    
    Ctj1 <- unique(Ctj[tempi.evento == tempo])
    Ctjx<-pp[tempi.evento == tempo,2:(num.sogg+1)]
    
    # for V_ind
    ris.V_num <- Ctjx*(1-Ctjx)
    ris.V_den <- (Ctj1*(1-Ctj1))
    
    # for V_w_ind
    ris.V_W <- Ctjx*(1-Ctjx) / (Ctj1*(1-Ctj1))
    
    # for DN
    rd <- ((Ctj1 - Ctjx) / Ctj1)^2
    ris.DN <- rd[(Ctjx<Ctj1) & (Ctj1>0)]

    # for DS
    rs <- ((Ctjx - Ctj1) / (1-Ctj1))^2
    ris.DS <- rs[(Ctjx>Ctj1) & (Ctj1<1)]
    
    return(c(ris.V_den, sum(ris.V_num)/num.sogg, 
             1-sum(ris.V_W)/num.sogg, 
             sqrt(sum(ris.DN)/sum((Ctjx<Ctj1) & (Ctj1>0))), 
             sqrt(sum(ris.DS)/sum((Ctjx>Ctj1) & (Ctj1<1))))) # sum gives i-loop in SAS
  }

  
  #function which calculates  M(t_(j)) for a given time t_(j) = tempo
  f.Mt <- function(tempo, tutti.tempi, 
                   Cti, tempi.evento, 
                   C_CE_ti,
                   Ctj, status, failcode, CEcode, num.sogg, 
                   Gtj, Gti)
  {
    Ctj1 <- unique(Ctj[tempi.evento == tempo])
    
    primo <- rep(Ctj1, num.sogg)
    primo[tutti.tempi <= tempo] <- 0
    pij <- (1-Ctj1-C_CE_ti)/(1-Cti-C_CE_ti)
    terzo <- (status==0) * (Ctj1*pij + (1-Ctj1)*(1 - pij))
    terzo[tutti.tempi > tempo] <- 0
    terzo[is.na(terzo)] <- 0
    secondo <- (1-Ctj1) * (status==failcode)
    secondo[tutti.tempi > tempo] <- 0
    quatro <- Ctj1 * (status==CEcode)
    quatro[tutti.tempi > tempo] <- 0
    ris <- primo + secondo + terzo + quatro
    return(sum(ris)/num.sogg) # sum gives i-loop
  }
  
  #function which calculates M(t_(j)|x) for a given time t_(j)
  f.Mt.fg <- function(tempo, tutti.tempi, 
                      Ctix, tempi.evento, 
                      C_CE_tix,
                      pp, status, failcode, CEcode, num.sogg, 
                      Gtj, Gti)
  {
    Ctjx<-pp[tempi.evento == tempo,2:(num.sogg+1)] # (1st column: tempi.eventi)
    
    primo <- Ctjx
    primo[tutti.tempi <= tempo] <- 0
    pij <- (1-Ctjx-C_CE_tix)/(1-Ctix-C_CE_tix)
    terzo <- (status==0) * (Ctjx*pij + (1-Ctjx)*(1 - pij))
    terzo[tutti.tempi > tempo] <- 0
    terzo[is.na(terzo)] <- 0
    secondo <- (1-Ctjx) * (status==failcode)
    secondo[tutti.tempi > tempo] <- 0
    quatro <- Ctjx * (status==CEcode)
    quatro[tutti.tempi > tempo] <- 0
    ris <- primo + secondo + terzo + quatro
    return(sum(ris)/num.sogg) # sum gives i-loop
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
  dat<-na.omit(cbind(time,status,x))
  time<-dat[,1]
  status<-dat[,2]
  x<-dat[,3]
  
  num.sogg <- length(time)
  
  cif.unc <- cuminc(time,status)
  
  # unconditional
  tempi.eventi <- sort(unique(time[status==failcode])) # different EoI times
  pos.cif <- apply(as.matrix(time), 1, f.assegna.surv, tempi.eventi) # position in tempi.eventi to read CIF at each individual's observed time
  unique.time.est<-timepoints(cif.unc,tempi.eventi)$est[failcode,]
  surv.tot.cif.unc <- unique.time.est[pos.cif] # each individual's uncond. CIF estimate at its observed time
  
  ind.censura <- as.numeric(!as.logical(status))
  cif.tj <- unique.time.est 
  numero.eventi <- as.vector(table(sort(time[status==failcode])))
  
  # conditional
  fgfit<-crr(time,status, x, failcode=failcode)
  pp <- predict(fgfit, x) # (1st column: tempi.eventi)
  cif.tot.fg <- pp[cbind(pos.cif,2:(num.sogg+1))]

  # unconditional for CE
  tempi.eventiCE <- sort(unique(time[status==CEcode])) # different CE times
  pos.cifCE <- apply(as.matrix(time), 1, f.assegna.surv, tempi.eventiCE) # position in tempi.eventiCE to read CIF_CE at each individual's observed time
  unique.time.estCE<-timepoints(cif.unc,tempi.eventiCE)$est[CEcode,]
  surv.tot.cifCE.unc <- unique.time.estCE[pos.cifCE] # each individual's uncond. CIF_CE estimate at its observed time
  
  # conditional for CE
  fgfitCE<-crr(time,status, x, failcode=CEcode)
  ppCE <- predict(fgfitCE, x) # (1st column: tempi.eventi)
  cifCE.tot.fg <- ppCE[cbind(pos.cifCE,2:(num.sogg+1))]

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
              surv.tot.cifCE.unc,
              cif.tj, status, failcode, CEcode, num.sogg,
              cens.tot.km, cens.tot.km[pos.cif]) # apply gives j-loop
  Mtx <- apply(as.matrix(tempi.eventi), 1, f.Mt.fg, time, 
               cif.tot.fg, tempi.eventi, 
               cifCE.tot.fg,
               pp, status, failcode, CEcode, num.sogg, 
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

  return(list(V = V, 
              V_ind = V_ind, Vw = Vw, Vw_ind = Vw_ind, 
              DN = DN1, DS = DS1))
}


  

### example

#library(cmprsk)

#cif.unc<-cuminc(abcsg$tt_dist, abcsg$stat_dist)
#  plot(cif.unc)
#fg.tusize<-crr(abcsg$tt_dist, abcsg$stat_dist, abcsg$log_TumorSizeN)
#  summary(fg.tusize)
#NecSuff_CR(abcsg$tt_dist, abcsg$stat_dist, abcsg$log_TumorSizeN)
