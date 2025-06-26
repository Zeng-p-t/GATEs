
rm(list = ls( ))
setwd('C:/Users/衍/OneDrive/桌面/Application')
library(locpol)
library(dplyr)
library(KernSmooth)   # used to select the optimal bandwidth
library(np)

# ----------------------  data preprocessing  ---------------------
load('data_longstudy.Rdata')    # 5034


days <- as.numeric(dat$FollowTime.1 - dat$BaseTime)
days = as.numeric(dat$FollowTime.1 - dat$BaseTime)
dat <- dat[days>30,]        # remove some samples

dat$FollowTime.1 <- NULL
dat$BaseTime <- NULL


# A <- dat$D    # treatment
dat$PASI.q936 <- ifelse(dat$PASI.q936 <= 0.25*dat$PASI.q155, 1, 0) # outcome, 75% PASI decrease
# Y <- dat$PASI.q936
colnames(dat)[c(14:15)] <- c('Y', 'A')  # outcome and treatment


# Design matrix
dat <- model.matrix(~ ., dat)
dat <- data.frame(dat)
dat$X.Intercept. <- NULL

# covariates
# X <- select(dat, -c('Y', 'A'))
# X <- as.matrix(X)

# ----------------------  data analysis  ---------------------

source('EstTau.R')
source('MatchY1Y0.R')
source('EstPSR.R')

############# step 1. fixed h ############################
n <- nrow(dat)
p <- 15
A <- dat$A                  # treatment
Y  <-  dat$Y                # outcome
X <- as.matrix(dat[, 1:p])  # covariates

Z <- dat$age
Zeval = seq(25, 65, by = 10)



Result <- data.frame(matrix(nrow = 1, ncol = 3))
Result[, 1] <- c('h')
names(Result) <- c('\040','match','match.bc') # ,'ipw','or','aipw', 'psr')

H <- data.frame(matrix(nrow = 500, ncol = 2))

set.seed(222)

for(i in 1:500){

    ########### generate data ###############
    n1 <- sum(A == 1)
    n0 <- sum(A == 0)
    index.treat <- which(A == 1)
    index.control <- which(A == 0)
    n1.s <- round(4*n1^(2/3))
    n0.s <- round(4*n0^(2/3))

    sub.i <- c(sample(index.treat, n1.s), sample(index.control, n0.s))
    dat.i <- dat[sub.i,]
    X.i <- as.matrix(dat.i[,1:p])
    A.i <- dat.i$A
    Y.i <- dat.i$Y
    Z.i <- dat.i$age

    ###### The bandwidth of the matching method #####
    res.match <- match_y1y0(X.i, A.i, Y.i, K = 5)   # one-to-five
    y1_y0.i <- res.match$Y1 - res.match$Y0

    h1.i <- robustdpill(x= Z.i, y = y1_y0.i)
    H[i, 1] <- h1.i



    ###### The bandwidth of the matching.bc method #####
    miu1_hat.i <- cbind(1,X.i)%*%as.matrix(lm(Y.i ~ X.i, subset = A.i==1)$coef)
    miu0_hat.i <- cbind(1,X.i)%*%as.matrix(lm(Y.i ~ X.i, subset = A.i==0)$coef)

    res.match <- match_y1y0_bc(X.i, A.i, Y.i,
                               miu1.hat = miu1_hat.i, miu0.hat = miu0_hat.i,
                               K = 5)   # one-to-five
    y1_y0.i <- res.match$Y1 - res.match$Y0

    h2.i <- robustdpill(x= Z.i, y = y1_y0.i)
    H[i, 2] <- h2.i

    ###### The bandwidth of the ipw method #####
    ###### The bandwidth of the aipw method #####
    ###### The bandwidth of the or method #####
    ###### The bandwidth of the psr method #####

    #cat('h1=',h1.i, ' h2=', h2.i,'\n')

    cat(i, '\r')
    i <- i + 1

}

Result[,2:3] <- apply(H, 2, mean, na.rm = T)
print(Result)


############# Step 2:  estimators variance #############

Variance <- data.frame(matrix(nrow = 6, ncol = 6))
Variance[, 1] <- c('match','match.bc','ipw','aipw','or','prs')
names(Variance) <- c('\040','25','35','45','55','65')

V1 <- data.frame(matrix(nrow = 500, ncol = 5))
V2 <- data.frame(matrix(nrow = 500, ncol = 5))
V3 <- data.frame(matrix(nrow = 500, ncol = 5))
V4 <- data.frame(matrix(nrow = 500, ncol = 5))
V5 <- data.frame(matrix(nrow = 500, ncol = 5))
V6 <- data.frame(matrix(nrow = 500, ncol = 5))

set.seed(666)

for(i in 1:500){

  ########### generate data ###############
  n1 <- sum(A == 1)
  n0 <- sum(A == 0)
  index.treat <- which(A == 1)
  index.control <- which(A == 0)
  n1.s <- round(4*n1^(2/3))
  n0.s <- round(4*n0^(2/3))

  sub.i <- c(sample(index.treat, n1.s), sample(index.control, n0.s))
  dat.i <- dat[sub.i,]
  X.i <- as.matrix(dat.i[,1:p])
  A.i <- dat.i$A
  Y.i <- dat.i$Y
  Z.i <- dat.i$age

  ## (1). matching-based estimator
  # obtaining Y1_hat - Y0_hat with matching
  res.match <- match_y1y0(X.i, A.i, Y.i, K = 5)   # one-to-five
  y1_y0.i <- res.match$Y1 - res.match$Y0
  # estimating GATE with local kernel regression
  h <- 5.68
  res <- EstTauCon(Y1_Y0 = y1_y0.i, Z.i, Zeval, h = h)
  est.match.i <- res$GATE
  V1[i, ] <- est.match.i

  ## (2). bias-corrected matching estimator
  # estimating outcome regression functions with linear models
  miu1_hat.i <- cbind(1,X.i)%*%as.matrix(lm(Y.i ~ X.i, subset = A.i==1)$coef)
  miu0_hat.i <- cbind(1,X.i)%*%as.matrix(lm(Y.i ~ X.i, subset = A.i==0)$coef)

  # obtaining bias-corrected Y1_hat - Y0_hat with matching
  res.match <- match_y1y0_bc(X.i, A.i, Y.i,
                             miu1.hat = miu1_hat.i, miu0.hat = miu0_hat.i,
                             K = 5)   # one-to-five
  y1_y0.i <- res.match$Y1 - res.match$Y0

  # estimating GATE with local kernel regression
  res <- EstTauCon(Y1_Y0 = y1_y0.i, Z.i, Zeval, h = h)
  est.match.bc.i <- res$GATE
  V2[i, ] <- est.match.bc.i

  # estimating propensity scores with logit regression
  ps_est.i <- fitted(glm(A.i~X.i, family=binomial("logit")))
  ps_est.i <- unname(ps_est.i)

  ##  IPW estimator
  psi1.i <- A.i*Y.i/ps_est.i
  psi0.i <- (1-A.i)*Y.i/(1-ps_est.i)
  psi.i <- psi1.i - psi0.i
  res <- EstTauCon(Y1_Y0 = psi.i, Z.i, Zeval, h = h)
  est.ipw.i <- res$GATE
  V3[i, ] <- est.ipw.i

  ##  OR estimator
  psi.i <- miu1_hat.i - miu0_hat.i
  res <- EstTauCon(Y1_Y0 = psi.i, Z.i, Zeval, h = h)
  est.or.i <- res$GATE
  V4[i, ] <- est.or.i


  ##  AIPW estimator
  psi1.i <- A.i*Y.i/ps_est.i - (A.i-ps_est.i)*miu1_hat.i/ps_est.i
  psi0.i <- (1-A.i)*Y.i/(1-ps_est.i)+(A.i-ps_est.i)*miu0_hat.i/(1-ps_est.i)
  psi.i <- psi1.i - psi0.i
  res <- EstTauCon(Y1_Y0 = psi.i, Z.i, Zeval, h = h)
  est.aipw.i <- res$GATE
  V5[i, ] <- est.aipw.i

  ##  PSR estimator
  h.i <- npscoefbw(formula = Y.i ~ A.i | Z.i + ps_est.i, nmulti = 1)$bw
  h1.i <- h.i[1]     # z
  h2.i <- h.i[2]      # ps
  res <- EstPSR(Y.i, A.i, ps_est.i, Z.i, Zeval, h1.i, h2.i, h3 = h.i)
  est.psr.i <- res$GATE
  V6[i, ] <- est.psr.i


  cat(i, '\r')
  i <- i + 1

}

Variance[1,2:6] <- apply(V1, 2, sd, na.rm = T)
Variance[2,2:6] <- apply(V2, 2, sd, na.rm = T)
Variance[3,2:6] <- apply(V3, 2, sd, na.rm = T)
Variance[4,2:6] <- apply(V4, 2, sd, na.rm = T)
Variance[5,2:6] <- apply(V5, 2, sd, na.rm = T)
Variance[6,2:6] <- apply(V6, 2, sd, na.rm = T)

print(Variance)
