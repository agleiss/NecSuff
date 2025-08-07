##################################################################
#
# NecSuff
# =======
#
# R-function for plotting the predictiveness curve and showing
# areas that correspond to DN and DS;
# 
# Gleiss, A. Visualizing a marker's degrees of necessity and of 
# sufficiency in the predictiveness curve, submitted (2024)
#
# Author:  Andreas Gleiss
# Version: 1.0
# Date:    17 Sept 2024
#
# Arguments:
# ==========
# 
# pred 			name of variable containing predictions
#           (output of predict())
#
##################################################################

NecSuff_PredCurve <- function(pred) {
  
  q.pred <- rank(pred, ties.method = "first") / length(pred)

  plt <- cbind(q.pred, pred)
  plts <- plt[order(q.pred),]

  P.D <- mean(pred)
  q0 <- mean(pred < P.D)

  par.orig <- par(pty = "s", # squared
                  omi = c(0, 0, 0, 0),
                  xaxs = "i", yaxs = "i") # 

  plot(NULL, xlim = c(0,1), ylim = c(0,1), 
       xlab = "Risk quantile q", ylab = "Risk R(q)")
  #abline(h = P.D, v = q0, col = "darkgrey")

  polygon(x = c(0, q0, q0, 0),
          y = c(0, 0, P.D, P.D), col = "lightgrey", border = NA)
  polygon(x = c(q0, 1, 1, q0),
          y = c(P.D, P.D, 1, 1), col = "lightgrey", border = NA)

  polygon(x = c(0, plts[,1], 1, 0, 0),
          y = c(0, plts[,2], P.D, P.D, 0), col = "darkgrey", border = NA)
  lines(plts)
  
  abline(v = c(0, 1), h = c(0, 1))
  
  par(par.orig)

}

#NecSuff(pred)
#NecSuff_PredCurve(pred)
