robustdpill <- function(x, y) {
  tryCatch(dpill(x,y),
           warning = function(w) {NaN},
           error = function(e){NaN})
}

# ------------  Estimating tau for continuous Z ------------
EstTauCon <- function(Y1_Y0, Z, Zeval, h){

  n <- length(Y1_Y0)

  ## select a bandwidth
  # h <-  robustdpill(x = W, y = beta1)  # select bandwidth by dpill
  # # h <-  h*n^(1/5)*n^(-2/7)
  # if (h == "NaN") {
  #   h <- npregbw(beta1 ~ W, regtype = 'll')$bw  #　select bandwidth via CV
  #   # h <- 3
  # }
  # print(h)

  tau <- locPolSmootherC(x=Z, y= Y1_Y0, xeval= Zeval,
                         bw= h, deg=0, kernel=gaussK)$beta0  #

  return( data.frame(Zeval = Zeval, GATE = tau))
}



# ------------  Estimating tau for discrete W ------------
# Note: W is a factor.
#
# EstTauDis <- function(beta1, W){
#
#   dat_tmp <- data.frame(beta1 = beta1, W =W)
#   by_w <- group_by(dat_tmp, W)
#   res <- summarize(by_w, tau = mean(beta1))
#
#   return( data.frame(res))
# }

