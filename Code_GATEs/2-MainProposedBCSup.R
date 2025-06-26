
rm(list = ls( ))
setwd("C:/Users/T14zpt/Desktop/Code")

library(np)
library(locpol)
library(dplyr)
# library(randomForest)

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

# ------------------------------------------- #
# case C1: c0 = 0.6 * n^(-1/5)
# case C2: c0 = 0.65 * n^(-1/5)
# case C3: c0 = 0.5 * n^(-1/5)
# 重复跑100次模拟, 平均窗宽约为c0 * n^(-1/5)
# 最后我们统一使用c0 = 0.6
# ------------------------------------------- #

set.seed(12)
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

    ## estimate outcome regression functions ##
    # mod.mu1 <- randomForest(formula = Y~.-A, data=dat,
    #                         subset = A == 1)
    # mod.mu0 <- randomForest(formula = Y~.-A, data=dat,
    #                         subset = A == 0)
    # miu1_hat <- unname(predict(mod.mu1, newdata = dat))
    # miu0_hat <- unname(predict(mod.mu0, newdata = dat))
    # mu0[sub.test] <- predict(mod.mu0, newdata = dat.test, type = 'prob')[,2]
    # X1 <- X[, 1]
    # X2 <- X[, 2]
    # X3 <- X[, 3]
    # XX <- cbind(X, X1^2, X1*X2, X3^3, X1*(2*X1+1)^2*(X1-1)^2, cos(3*X1) * exp(X1) * log(X1+2))
    # miu1_hat <- cbind(1,XX)%*%as.matrix(lm(Y ~ XX, subset = A==1)$coef)
    # miu0_hat <- cbind(1,XX)%*%as.matrix(lm(Y ~ XX, subset = A==0)$coef)
    miu1_hat <- cbind(1,X)%*%as.matrix(lm(Y ~ X, subset = A==1)$coef)
    miu0_hat <- cbind(1,X)%*%as.matrix(lm(Y ~ X, subset = A==0)$coef)

    ## obtaining Y1_hat - Y0_hat with matching
    res.match <- match_y1y0_bc(X, A, Y,
                                miu1.hat = miu1_hat, miu0.hat = miu0_hat,
                                K = 5)   # one-to-five
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
    n1.s <- round(3.3*n1^(2/3))
    n0.s <- round(3.3*n0^(2/3))

    GATE.B <- matrix(nrow = B, ncol = 5)
    for(k in 1:B){

      sub.k <- c(sample(index.treat, n1.s), sample(index.control, n0.s))
      # sub.k <- sample(n, replace = T)
      dat.k <- dat[sub.k, ]
      X.k <- as.matrix(dat.k[,1:p])
      A.k <- dat.k$A
      Y.k <- dat.k$Y
      Z.k <- dat.k$X1

      # mod.mu1 <- randomForest(formula = Y~.-A, data=dat.k,
      #                         subset = A == 1)
      # mod.mu0 <- randomForest(formula = Y~.-A, data=dat.k,
      #                         subset = A == 0)
      # miu1_hat.k <- unname(predict(mod.mu1, newdata = dat.k))
      # miu0_hat.k <- unname(predict(mod.mu0, newdata = dat.k))
        
      # X1 <- X.k[, 1]
      # X2 <- X.k[, 2]
      # X3 <- X.k[, 3]
      # XX.k <- cbind(X.k, X1^2, X1*X2, X3^3, X1*(2*X1+1)^2*(X1-1)^2, cos(3*X1) * exp(X1) * log(X1+2))
      # miu1_hat.k <- cbind(1,XX.k)%*%as.matrix(lm(Y.k ~ XX.k, subset = A.k==1)$coef)
      # miu0_hat.k <- cbind(1,XX.k)%*%as.matrix(lm(Y.k ~ XX.k, subset = A.k==0)$coef)

      miu1_hat.k <- cbind(1,X.k)%*%as.matrix(lm(Y.k ~ X.k, subset = A.k==1)$coef)
      miu0_hat.k <- cbind(1,X.k)%*%as.matrix(lm(Y.k ~ X.k, subset = A.k==0)$coef)
        
      res.match.k <- match_y1y0_bc(X.k, A.k, Y.k,
                                  miu1.hat = miu1_hat.k, miu0.hat = miu0_hat.k,
                                  K = 5)
      y1_y0.k <- res.match.k$Y1 - res.match.k$Y0

      GATE.B[k, ] <-  locPolSmootherC(x=Z.k, y= y1_y0.k, xeval = Zeval,
                                        bw= h, deg=0, kernel=gaussK)$beta0
      cat(k, '\r')
    }
    ese <- apply(GATE.B, 2, sd, na.rm = T)
    ese <- apply(GATE.B, 2, function(x){
      x <- sort(x); x <- sort(x)[5:96]; return(sd(x))
    })
    ESE[i, ] <- ese

    tau.true <- TrueTau[[ss]](Zeval)
    count <- (tau.est - 1.96*ese <= tau.true) & (tau.est + 1.96*ese >= tau.true)
    COUNT95[i, ] <- count

    cat(i, '\r')
    i <- i + 1
  }
  save(PARA, file = paste0('res_proposed/PARA_ss',ss, '_n=',n,'_bc.Rdata'))
  #save(ESE, file = paste0('result/PARA_ss',ss, '_n=',n,'_c=',c0,'_match.Rdata'))
  #save(COUNT95, file = paste0('result/PARA_ss',ss, '_n=',n,'_c=',c0,'_match.Rdata'))

  # report some results for z = -0.4, -0.2, 0.0, 0.2, 0.4
  RES <- data.frame(matrix(nrow = 5, ncol = 5))
  names(RES) <- c('z', 'match.bias','sd','ese','Cp95')

  RES$z <- Zeval
  RES[,2] <- apply(PARA, 2, mean, na.rm = T) - tau.true
  RES[,3] <- apply(PARA, 2, sd, na.rm = T)
  RES[,4] <- apply(ESE, 2, mean, na.rm = T)
  RES[,5] <- apply(COUNT95, 2, mean, na.rm = T)


  write.csv(RES, file = paste0('res_proposed/bc_ss',ss, '_n=',n,'.csv'),
            row.names = F)
}

