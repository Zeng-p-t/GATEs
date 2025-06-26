
rm(list = ls( ))
setwd("C:/Users/T14zpt/Desktop/Code")

library(np)
library(locpol)
library(dplyr)

# Matching-based estimator.
source('EstTau.R')
source('MatchY1Y0.R')
source('EstPSR.R')

# Functions for generating data
source('GenData.R')
GenData <- list(gen.dat1, gen.dat2, gen.dat3,
                gen.dat4, gen.dat5, gen.dat6,
                gen.dat7, gen.dat8, gen.dat9,
                gen.dat10, gen.dat11, gen.dat12)
TrueMu1 <- list(true.tau1, true.tau2, true.tau3,
                true.tau4, true.tau5, true.tau6,
                true.tau7, true.tau8, true.tau9,
                true.tau10, true.tau11, true.tau12)


set.seed(12)
for(ss in 1:12){

  LOOP = 1000
  PARA.match <- matrix(nrow = LOOP, ncol = 99)
  PARA.match.bc <- matrix(nrow = LOOP, ncol = 99)
  PARA.ipw <- matrix(nrow = LOOP, ncol = 99)
  PARA.or <- matrix(nrow = LOOP, ncol = 99)
  PARA.aipw <- matrix(nrow = LOOP, ncol = 99)
  PARA.psr <- matrix(nrow = LOOP, ncol = 99)

  n <- 2000
  p <- 3
  cat('====== ss-',ss, '  n =', n,  '  ===== \n')
  j <- 1

  while(j <= LOOP){  # LOOP

    ############# Step 1: generate data  ##############
    dat <- GenData[[ss]](n = n)

    A <- dat$A                  # treatment
    Y <- dat$Y                # outcome
    X <- as.matrix(dat[, 1:p])  # covariates

    h <- 0.06          # bandwidth for local constant regression
    Z <- dat$X1
    Zeval = seq(-0.49, 0.49, len = 99)

    ############# Step 2: proposed estimators #############

    ## (1). matching-based estimator
    # obtaining Y1_hat - Y0_hat with matching
    res.match <- match_y1y0(X, A, Y, K = 5)   # one-to-five
    y1_y0 <- res.match$Y1 - res.match$Y0
    # estimating GATE with local kernel regression
    res <- EstTauCon(Y1_Y0 = y1_y0, Z, Zeval, h = h)
    PARA.match[j, ] <- res$GATE

    ## (2). bias-corrected matching estimator
    # estimating outcome regression functions with linear models
    miu1_hat <- cbind(1,X)%*%as.matrix(lm(Y ~ X, subset = A==1)$coef)
    miu0_hat <- cbind(1,X)%*%as.matrix(lm(Y ~ X, subset = A==0)$coef)

    # obtaining bias-corrected Y1_hat - Y0_hat with matching
    res.match <- match_y1y0_bc(X, A, Y, miu1_hat, miu0_hat, K = 5)   # one-to-five
    y1_y0 <- res.match$Y1 - res.match$Y0
    # estimating GATE with local kernel regression
    res <- EstTauCon(Y1_Y0 = y1_y0, Z, Zeval, h = h)
    PARA.match.bc[j, ] <- res$GATE


    ############# Step 3: competing estimators #############

    # estimating propensity scores with logit regression
    ps_est <- fitted(glm(A~X, family=binomial("logit")))
    ps_est <- unname(ps_est)

    ## (1) IPW estimator
    psi1 <- A*Y/ps_est
    psi0 <- (1-A)*Y/(1-ps_est)
    psi <- psi1 - psi0
    res <- EstTauCon(Y1_Y0 = psi, Z, Zeval, h = h)
    PARA.ipw[j, ] <- res$GATE

    ## (2) OR estimator
    psi <- miu1_hat - miu0_hat
    res <- EstTauCon(Y1_Y0 = psi, Z, Zeval, h = h)
    PARA.or[j, ] <- res$GATE

    ## (3) AIPW estimator
    psi1 <- A*Y/ps_est - (A-ps_est)*miu1_hat/ps_est
    psi0 <- (1-A)*Y/(1-ps_est)+(A-ps_est)*miu0_hat/(1-ps_est)
    psi <- psi1 - psi0
    res <- EstTauCon(Y1_Y0 = psi, Z, Zeval, h = h)
    PARA.aipw[j, ] <- res$GATE

    ## (4) PSR estimator
    h1 <- 0.1
    h2 <- 0.1
    # h1 = 0.1 * n^{-1/6}
    # h2 = 0.1 * n^{-1/6}
    # h <- npscoefbw(formula = Y ~ A | Z + ps_est, nmulti = 1)$bw
    # h <- ifelse(h <= 0.1 * n^(-1/6), 0.1 * n^(-1/6), h)
    # h <- ifelse(h >=  n^(-1/6),   n^(-1/6), h)
    # h1 <- h[1]      # w
    # h2 <- h[2]      # ps
    res <- EstPSR(Y, A, ps_est, Z, Zeval, h1, h2, h3 = h)
    PARA.psr[j, ] <- res$GATE

    cat(j, '\r')
    j <- j + 1
  }
  save(PARA.match, file = paste0('result/PARA_ss',ss, '_n=',n,'_match.Rdata'))
  save(PARA.match.bc, file = paste0('result/PARA_ss',ss, '_n=',n,'_match.bc.Rdata'))
  save(PARA.ipw, file = paste0('result/PARA_ss',ss, '_n=',n,'_ips.Rdata'))
  save(PARA.or, file = paste0('result/PARA_ss',ss, '_n=',n,'_or.Rdata'))
  save(PARA.aipw, file = paste0('result/PARA_ss',ss, '_n=',n,'_aipw.Rdata'))
  save(PARA.psr, file = paste0('result/PARA_ss',ss, '_n=',n,'_psr.Rdata'))

  # report some results for z = -0.4, -0.2, 0.0, 0.2, 0.4 
  RES <- data.frame(matrix(nrow = 5, ncol = 2*6+1))
  names(RES) <- c('z', 'match.bias','sd', 'match.bc.bias', 'sd',
                    'ipw.bias', 'sd', 'or.bias', 'sd',
                    'aipw.bias', 'sd', 'psr.bias', 'sd')
  RES$z <- c(-0.4, -0.2, 0.0, 0.2, 0.4)
  index <- c(10,30,50,70,90)
  true.tau <- TrueMu1[[ss]](c(-0.4, -0.2, 0.0, 0.2, 0.4))
  RES[,2] <- apply(PARA.match[,index], 2, mean, na.rm = T) - true.tau
  RES[,3] <- apply(PARA.match[,index], 2, sd, na.rm = T)
  RES[,4] <- apply(PARA.match.bc[,index], 2, mean, na.rm = T) - true.tau
  RES[,5] <- apply(PARA.match.bc[,index], 2, sd, na.rm = T)
  RES[,6] <- apply(PARA.ipw[,index], 2, mean, na.rm = T) - true.tau
  RES[,7] <- apply(PARA.ipw[,index], 2, sd, na.rm = T)
  RES[,8] <- apply(PARA.or[,index], 2, mean, na.rm = T) - true.tau
  RES[,9] <- apply(PARA.or[,index], 2, sd, na.rm = T)
  RES[,10] <- apply(PARA.aipw[,index], 2, mean, na.rm = T) - true.tau
  RES[,11] <- apply(PARA.aipw[,index], 2, sd, na.rm = T)
  RES[,12] <- apply(PARA.psr[,index], 2, mean, na.rm = T) - true.tau
  RES[,13] <- apply(PARA.psr[,index], 2, sd, na.rm = T)

  write.csv(RES, file = paste0('result/result_ss',ss, '_n=',n,'.csv'),
             row.names = F)


}

