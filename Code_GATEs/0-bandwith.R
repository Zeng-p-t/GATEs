
rm(list = ls( ))
setwd("C:/Users/T14zpt/Desktop/Code")

library(KernSmooth)
library(dplyr)

robustdpill <- function(x, y) {
  tryCatch(dpill(x,y),
           warning = function(w) {NaN},
           error = function(e){NaN})
}


# Functions for estimating the conditional treatment effects.
# source('EstTau.R')
source('MatchY1Y0.R')

# Functions for generating data
source('GenData.R')
GenData <- list(gen.dat1, gen.dat2, gen.dat3,
                gen.dat4, gen.dat5, gen.dat6,
                gen.dat7, gen.dat8, gen.dat9)
TrueTau <- list(true.tau1, true.tau2, true.tau3,
                true.tau4, true.tau5, true.tau6,
                true.tau7, true.tau8, true.tau9)


H.match <- matrix(nrow = 100, ncol = 9)
H.match.bc <- matrix(nrow = 100, ncol = 9)  


set.seed(12)
for(ss in 1:9){
  
    LOOP <- 100
    
    
    n <- 2000
    p <- 3
    cat('====== ss-',ss, '  n =', n, '  ===== \n')
    
    for(i in 1:LOOP){  # LOOP
      
      ############# Step 1: generate data  ##############
      dat <- GenData[[ss]](n = n)
      
      ############# Step 2: conduct matching #############
      A <- dat$A   # treatment
      Y <- dat$Y
      X <- as.matrix(dat[, 1:p])
      Z <- dat$X1
      
      ## obtaining Y1_hat - Y0_hat with matching
      res.match <- match_y1y0(X, A, Y, K = 5)   # one-to-five
      y1_y0 <- res.match$Y1 - res.match$Y0
      
      h1 <- robustdpill(x= Z, y = y1_y0)
      H.match[i, ss] <- h1
      
      
      
      ## obtaining Y1_hat - Y0_hat with bias-corrected matching
      miu1_hat <- cbind(1,X)%*%as.matrix(lm(Y ~ X, subset = A==1)$coef)
      miu0_hat <- cbind(1,X)%*%as.matrix(lm(Y ~ X, subset = A==0)$coef)
      
      res.match <- match_y1y0_bc(X, A, Y,
                                 miu1.hat = miu1_hat, miu0.hat = miu0_hat,
                                 K = 5)   # one-to-five
      y1_y0 <- res.match$Y1 - res.match$Y0
  
      h2 <- robustdpill(x= Z, y = y1_y0) 
      H.match.bc[i, ss] <- h2
      
      cat('h1=',h1, ' h2=', h2, '\n')
      
      cat(i, '\r')
      i <- i + 1
    }
    
}

apply(H.match, 2, mean, na.rm = T)
apply(H.match.bc, 2, mean, na.rm = T)


