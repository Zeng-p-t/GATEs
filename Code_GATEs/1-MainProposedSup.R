
rm(list = ls( ))
setwd("C:/Users/T14zpt/Desktop/Code")

library(np)
library(locpol)
library(dplyr)

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


set.seed(123)
for(ss in 1:9){
  LOOP = 1000
  PARA <- matrix(nrow = LOOP, ncol = 5)
  ESE <- matrix(nrow = LOOP, ncol = 5)
  COUNT95 <- matrix(nrow = LOOP, ncol = 5)

  n <- 2000
  p <- 3
  cat('====== ss-',ss, '  n =', n, '  ===== \n')
  i <- 1

  while(i <= LOOP){  # LOOP

    ############# Step 1: generate data  ##############
    dat <- GenData[[ss]](n = n)

    ############# Step 2: conduct matching #############
    A <- dat$A   # treatment
    Y <- dat$Y
    X <- as.matrix(dat[, 1:p])

    ## obtaining Y1_hat - Y0_hat with matching
    res.match <- match_y1y0(X, A, Y, K = 5)   # one-to-five
    y1_y0 <- res.match$Y1 - res.match$Y0

    #### Step 3: estimate the heterogeneous treatment effects #####
    Z <- dat$X1
    Zeval = c(-0.4, -0.2, 0, 0.2, 0.4)

    h <- 0.06

    tau.est <- locPolSmootherC(x=Z, y= y1_y0, xeval = Zeval,
                             bw= h, deg=0, kernel=gaussK)$beta0
    PARA[i, ] <- tau.est

    #### estimating ESE ####
    B <- 100
    n1 <- sum(A == 1)
    n0 <- sum(A == 0)
    index.treat <- which(A == 1)
    index.control <- which(A == 0)
    n1.s <- round(4.15*n1^(2/3))
    n0.s <- round(4.15*n0^(2/3))

    GATE.B <- matrix(nrow = B, ncol = 5)
    for(k in 1:B){

      sub.k <- c(sample(index.treat, n1.s), sample(index.control, n0.s))
      dat.k <- dat[sub.k,]
      X.k <- as.matrix(dat.k[,1:p])
      A.k <- dat.k$A
      Y.k <- dat.k$Y
      Z.k <- dat.k$X1

      res.match.k <- match_y1y0(X.k, A.k, Y.k, K = 5)
      y1_y0.k <- res.match.k$Y1 - res.match.k$Y0

      GATE.B[k, ] <-  locPolSmootherC(x=Z.k, y= y1_y0.k, xeval = Zeval,
                                      bw= h, deg=0, kernel=gaussK)$beta0
      cat(k, '\r')
    }
    ese <- apply(GATE.B, 2, sd, na.rm = T)
    ESE[i, ] <- ese


    tau.true <- TrueTau[[ss]](Zeval)
    count <- (tau.est - 1.96*ese <= tau.true) & (tau.est + 1.96*ese >= tau.true)
    COUNT95[i, ] <- count

    cat(i, '\r')
    i <- i + 1
  }
  save(PARA, file = paste0('res_proposed/PARA_ss',ss, '_n=',n,'_proposed.Rdata'))
  save(ESE, file = paste0('result/PARA_ss',ss, '_n=',n,'_c=',c0,'_match.Rdata'))
  save(COUNT95, file = paste0('result/PARA_ss',ss, '_n=',n,'_c=',c0,'_match.Rdata'))

  # report some results for z = -0.4, -0.2, 0.0, 0.2, 0.4
  RES <- data.frame(matrix(nrow = 5, ncol = 5))
  names(RES) <- c('z', 'match.bias','sd','ese','Cp95')

  RES$z <- Zeval
  RES[,2] <- apply(PARA, 2, mean, na.rm = T) - tau.true
  RES[,3] <- apply(PARA, 2, sd, na.rm = T)
  RES[,4] <- apply(ESE, 2, mean, na.rm = T)
  RES[,5] <- apply(COUNT95, 2, mean, na.rm = T)

  print(RES)

  write.csv(RES, file = paste0('res_proposed/ss',ss, '_n=',n,'_sup.csv'),
              row.names = F)

    
}

