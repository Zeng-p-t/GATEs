robustdpill <- function(x, y) {
  tryCatch(dpill(x,y),
           warning = function(w) {NaN},
           error = function(e){NaN})
}


lm_C <- function (x, y, w, tol = 1e-07) {
  
  zero.weights <- any(w == 0)
  if (zero.weights) {
    ok <- w != 0
    w <- w[ok]
    x <- x[ok, , drop = FALSE]
    y <-  y[ok]
  }
  p <- ncol(x)
  wts <- sqrt(w)
  z <- .Call(stats:::C_Cdqrls, x * wts, y * wts, tol, FALSE)
  
  coef <- z$coefficients
  return(coef)
}



EstPSR <- function(Y, A, Pscore, W, Weval, h1, h2, h3){
  
  # 1. Estimate g1(e,w) and g2(e,w)
  n <- length(A)

  Z <- data.frame( t(cbind(W, Pscore)) )
  Beta <- sapply( Z, function(z){
    tmp1 <- (W - z[1])/h1; tmp2 <-  (Pscore - z[2])/h2
    
    w <- dnorm(tmp1) * dnorm(tmp2) / (h1 * h2)      # gaussian kernel
    w <- ifelse(is.na(w) , 0, w)
    xx <- cbind(A, 1) #D*tmp1, tmp1,D*tmp2, tmp2)
    coef <- lm_C(xx, Y, w)[1:2]
    return( coef )
  })
  Beta <- as.data.frame(t(Beta))
  names(Beta) <- c('beta1', 'beta2')
  
  g1_ll <- Beta$beta1
  g1_lc <- Beta$beta2
  xi <- Y -  A*g1_ll - g1_lc       # which is equal to Y - mod_ll$mean
  
  # 2. Estimate E[g1(e,w)|w]
  # 2.1 select a bandwidth
  # h3 <-  robustdpill(x = W, y = g1_ll)       # select bandwidth by dpill
  # if (h3 == "NaN") {h3 <- 0.1}

  tau <- locPolSmootherC(x=W, y= g1_ll, xeval= Weval,
                         bw= h3, deg=0, kernel=gaussK)$beta0

  return(data.frame(Zeval = Weval, GATE = tau))
}


