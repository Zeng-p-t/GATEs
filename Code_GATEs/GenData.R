# -------------  Cases (1) - (3)  -------------------
# both PS and OR are misspecified, 
# Propensity scores are not close to 0 or 1. 

gen.dat1 <- function(n){
  
  # X = (X1, X2, X3)
  X1 <- runif(n, -0.5,0.5) # rnorm(n)   
  X2 <- sample(c(0,1,2), n, T)# rnorm(n)
  X3 <- rnorm(n)
    
  # A 
  alpha <-  c(1/2, 1/4, -1/8)  # rep(0.5, p)
  # alpha <-  c(1, 1, -1, -1) * 0.25
  temp <-  exp( cbind(X1^2, X2^2, X3^2)  %*% alpha)
  A <-  ifelse(runif(n) <=   temp / (1+temp), 1, 0)  
  # hist(temp / (1+temp), nclass = 50)
  
  # Y0, Y1
  Y0 <-    X2 + X1*X2 + 0.5*(X3^3 + X3) + rnorm(n)
  Y1 <- A * (2*X1^2) + X2 + X1*X2 + 0.5*(X3^3 + X3) + rnorm(n)    # Yt
  Y <- A*Y1 + (1-A)*Y0 
  
  dat <- data.frame(X1 = X1, X2 = X2, X3 = X3, A = A, Y = Y)
  return(dat)
}




gen.dat2 <- function(n){
  
  # X = (X1, X2, X3)
  X1 <- runif(n, -0.5,0.5) # rnorm(n)   
  X2 <- sample(c(0,1,2), n, T)# rnorm(n)
  X3 <- rnorm(n)
  
  # A 
  # alpha <-  c(1, 1/2, -1/4, -1/8)
  alpha <-  c(1/2, 1/4, -1/8)  # rep(0.5, p)
  # alpha <-  c(1, 1, -1, -1) * 0.5
  temp <-  exp( cbind(X1^2, X2^2, X3^2)  %*% alpha)
  # alpha <-  c(1, 1, -1) * 0.25  # rep(0.5, p)
  # temp <-  exp( cbind(X1^2, X2^2, X3^2)  %*% alpha)
  A <-  ifelse(runif(n) <=   temp / (1+temp), 1, 0)  
  
  # Y0, Y1
  Y0 <-  X2 + X1*X2 + 0.5*(X3^3 + X3) + rnorm(n)
  Y1 <- A * (X1*(2*X1+1)^2*(X1-1)^2) + X2 + X1*X2 + 0.5*(X3^3 + X3) + rnorm(n)    # Yt
  Y <- A*Y1 + (1-A)*Y0 
  
  dat <- data.frame(X1 = X1, X2 = X2, X3 = X3, A = A, Y = Y)
  return(dat)
}




gen.dat3 <- function(n){
  
  # X = (X1, X2, X3)
  X1 <- runif(n, -0.5,0.5) # rnorm(n)   
  X2 <- sample(c(0,1,2), n, T)# rnorm(n)
  X3 <- rnorm(n)
  
  # A 
  alpha <-  c(1/2, 1/4, -1/8)  # rep(0.5, p)
  # alpha <-  c(1, 1, -1, -1) * 0.25
  temp <-  exp( cbind(X1^2, X2^2, X3^2)  %*% alpha)
  # alpha <-  c(1, 1, -1) * 0.25  # rep(0.5, p)
  # temp <-  exp( cbind(X1^2, X2^2, X3^2)  %*% alpha)
  A <-  ifelse(runif(n) <=   temp / (1+temp), 1, 0)  
  
  # Y0, Y1
  Y0 <-   X2 + X1*X2 + 0.5*(X3^3 + X3) + rnorm(n)
  Y1 <- A * (cos(3*X1) * exp(X1) * log(X1+2)) + X2 + X1*X2 + 0.5*(X3^3 + X3) + rnorm(n)   
  Y <- A*Y1 + (1-A)*Y0 
  
  dat <- data.frame(X1 = X1, X2 = X2, X3 = X3, A = A, Y = Y)
  return(dat)
}


# -------------  Cases (4) - (6)  -------------------
# both PS and OR are misspecified, 
# there are some propensity scores that close to 0 or 1. 

gen.dat4 <- function(n){
  
  # X = (X1, X2, X3)
  X1 <- runif(n, -0.5,0.5) # rnorm(n)   
  X2 <- sample(c(0,1,2), n, T)# rnorm(n)
  X3 <- rnorm(n)
  
  # A 
  alpha <-  c(1/2, 1/4, -1/8)*2  # rep(0.5, p)
  # alpha <-  c(1, 1, -1, -1) * 0.25
  temp <-  exp( cbind(8*X1^2, X2^2, 5*X3^2)  %*% alpha)
  A <-  ifelse(runif(n) <=   temp / (1+temp), 1, 0)  
  # hist(temp / (1+temp), nclass = 50)
  
  # Y0, Y1
  Y0 <-    X2 + X1*X2 + 0.5*(X3^3 + X3) + rnorm(n)
  Y1 <- A * (2*X1^2) + X2 + X1*X2 + 0.5*(X3^3 + X3) + rnorm(n)    # Yt
  Y <- A*Y1 + (1-A)*Y0 
  
  dat <- data.frame(X1 = X1, X2 = X2, X3 = X3, A = A, Y = Y)
  return(dat)
}




gen.dat5 <- function(n){
  
  # X = (X1, X2, X3)
  X1 <- runif(n, -0.5,0.5) # rnorm(n)   
  X2 <- sample(c(0,1,2), n, T)# rnorm(n)
  X3 <- rnorm(n)
  
  # A 
  # alpha <-  c(1, 1/2, -1/4, -1/8)
  alpha <-  c(1/2, 1/4, -1/8)*2  # rep(0.5, p)
  # alpha <-  c(1, 1, -1, -1) * 0.5
  temp <-  exp( cbind(8*X1^2, X2^2, 5*X3^2)  %*% alpha)
  # alpha <-  c(1, 1, -1) * 0.25  # rep(0.5, p)
  # temp <-  exp( cbind(X1^2, X2^2, X3^2)  %*% alpha)
  A <-  ifelse(runif(n) <=   temp / (1+temp), 1, 0)  
  
  # Y0, Y1
  Y0 <-  X2 + X1*X2 + 0.5*(X3^3 + X3) + rnorm(n)
  Y1 <- A * (X1*(2*X1+1)^2*(X1-1)^2) + X2 + X1*X2 + 0.5*(X3^3 + X3) + rnorm(n)    # Yt
  Y <- A*Y1 + (1-A)*Y0 
  
  dat <- data.frame(X1 = X1, X2 = X2, X3 = X3, A = A, Y = Y)
  return(dat)
}


gen.dat6 <- function(n){
  
  # X = (X1, X2, X3)
  X1 <- runif(n, -0.5,0.5) # rnorm(n)   
  X2 <- sample(c(0,1,2), n, T)# rnorm(n)
  X3 <- rnorm(n)
  
  # A 
  alpha <-  c(1/2, 1/4, -1/8)*2  # rep(0.5, p)
  # alpha <-  c(1, 1, -1, -1) * 0.25
  temp <-  exp( cbind(8*X1^2, X2^2, 5*X3^2)  %*% alpha)
  # alpha <-  c(1, 1, -1) * 0.25  # rep(0.5, p)
  # temp <-  exp( cbind(X1^2, X2^2, X3^2)  %*% alpha)
  A <-  ifelse(runif(n) <=   temp / (1+temp), 1, 0)  
  
  # Y0, Y1
  Y0 <-   X2 + X1*X2 + 0.5*(X3^3 + X3) + rnorm(n)
  Y1 <- A * (cos(3*X1) * exp(X1) * log(X1+2)) + X2 + X1*X2 + 0.5*(X3^3 + X3) + rnorm(n)   
  Y <- A*Y1 + (1-A)*Y0 
  
  dat <- data.frame(X1 = X1, X2 = X2, X3 = X3, A = A, Y = Y)
  return(dat)
}



# -------------  Cases (7) - (9)  -------------------
# PS is correctly specified, but OR is misspecified, 
# there are some propensity scores that close to 0 or 1. 

gen.dat7 <- function(n){
  
  # X = (X1, X2, X3)
  X1 <- runif(n, -0.5,0.5) # rnorm(n)   
  X2 <- sample(c(0,1,2), n, T)# rnorm(n)
  X3 <- rnorm(n)
  
  # A 
  alpha <-  c(3/2, 3/2, 3/2)  # rep(0.5, p)
  # alpha <-  c(1, 1, -1, -1) * 0.25
  temp <-  exp( cbind(X1, X2, X3)  %*% alpha)
  A <-  ifelse(runif(n) <=   temp / (1+temp), 1, 0)  
  # hist(temp / (1+temp), nclass = 50)
  
  # Y0, Y1
  Y0 <-    X2 + X1*X2 + 0.5*(X3^3 + X3) + rnorm(n)
  Y1 <- A * (2*X1^2) + X2 + X1*X2 + 0.5*(X3^3 + X3) + rnorm(n)    # Yt
  Y <- A*Y1 + (1-A)*Y0 
  
  dat <- data.frame(X1 = X1, X2 = X2, X3 = X3, A = A, Y = Y)
  return(dat)
}




gen.dat8 <- function(n){
  
  # X = (X1, X2, X3)
  X1 <- runif(n, -0.5,0.5) # rnorm(n)   
  X2 <- sample(c(0,1,2), n, T)# rnorm(n)
  X3 <- rnorm(n)
  
  # A 
  # alpha <-  c(1, 1/2, -1/4, -1/8)
  alpha <-  c(3/2, 3/2, 3/2)  # rep(0.5, p)
  # alpha <-  c(1, 1, -1, -1) * 0.5
  temp <-  exp( cbind(X1, X2, X3)  %*% alpha)
  # alpha <-  c(1, 1, -1) * 0.25  # rep(0.5, p)
  # temp <-  exp( cbind(X1^2, X2^2, X3^2)  %*% alpha)
  A <-  ifelse(runif(n) <=   temp / (1+temp), 1, 0)  
  
  # Y0, Y1
  Y0 <-  X2 + X1*X2 + 0.5*(X3^3 + X3) + rnorm(n)
  Y1 <- A * (X1*(2*X1+1)^2*(X1-1)^2) + X2 + X1*X2 + 0.5*(X3^3 + X3) + rnorm(n)    # Yt
  Y <- A*Y1 + (1-A)*Y0 
  
  dat <- data.frame(X1 = X1, X2 = X2, X3 = X3, A = A, Y = Y)
  return(dat)
}




gen.dat9 <- function(n){
  
  # X = (X1, X2, X3)
  X1 <- runif(n, -0.5,0.5) # rnorm(n)   
  X2 <- sample(c(0,1,2), n, T)# rnorm(n)
  X3 <- rnorm(n)
  
  # A 
  alpha <-  c(3/2, 3/2, 3/2)  # rep(0.5, p)
  # alpha <-  c(1, 1, -1, -1) * 0.25
  temp <-  exp( cbind(X1, X2, X3)  %*% alpha)
  # alpha <-  c(1, 1, -1) * 0.25  # rep(0.5, p)
  # temp <-  exp( cbind(X1^2, X2^2, X3^2)  %*% alpha)
  A <-  ifelse(runif(n) <=   temp / (1+temp), 1, 0)  
  
  # Y0, Y1
  Y0 <-   X2 + X1*X2 + 0.5*(X3^3 + X3) + rnorm(n)
  Y1 <- A * (cos(3*X1) * exp(X1) * log(X1+2)) + X2 + X1*X2 + 0.5*(X3^3 + X3) + rnorm(n)   
  Y <- A*Y1 + (1-A)*Y0 
  
  dat <- data.frame(X1 = X1, X2 = X2, X3 = X3, A = A, Y = Y)
  return(dat)
}



# -------------  Cases (10) - (12)  -------------------
# the true propensity scores can take values very close to 0.

gen.dat10 <- function(n){
  
  # X = (X1, X2, X3)
  X1 <- runif(n, -0.5,0.5) # rnorm(n)   
  X2 <- sample(c(0,1,2), n, T)# rnorm(n)
  X3 <- rnorm(n)
  
  # A 
  alpha <-  c(1/2, 1/4, -1/8)*3  # rep(0.5, p)
  # alpha <-  c(1, 1, -1, -1) * 0.25
  temp <-  exp( cbind(8*X1^2, X2^2, 5*X3^2)  %*% alpha)
  A <-  ifelse(runif(n) <=   temp / (1+temp), 1, 0)  
  # hist(temp / (1+temp), nclass = 50)
  
  # Y0, Y1
  Y0 <-    X2 + X1*X2 + 0.5*(X3^3 + X3) + rnorm(n)
  Y1 <- A * (2*X1^2) + X2 + X1*X2 + 0.5*(X3^3 + X3) + rnorm(n)    # Yt
  Y <- A*Y1 + (1-A)*Y0 
  
  dat <- data.frame(X1 = X1, X2 = X2, X3 = X3, A = A, Y = Y)
  return(dat)
}




gen.dat11 <- function(n){
  
  # X = (X1, X2, X3)
  X1 <- runif(n, -0.5,0.5) # rnorm(n)   
  X2 <- sample(c(0,1,2), n, T)# rnorm(n)
  X3 <- rnorm(n)
  
  # A 
  # alpha <-  c(1, 1/2, -1/4, -1/8)
  alpha <-  c(1/2, 1/4, -1/8)*3  # rep(0.5, p)
  # alpha <-  c(1, 1, -1, -1) * 0.5
  temp <-  exp( cbind(8*X1^2, X2^2, 5*X3^2)  %*% alpha)
  # alpha <-  c(1, 1, -1) * 0.25  # rep(0.5, p)
  # temp <-  exp( cbind(X1^2, X2^2, X3^2)  %*% alpha)
  A <-  ifelse(runif(n) <=   temp / (1+temp), 1, 0)  
  
  # Y0, Y1
  Y0 <-  X2 + X1*X2 + 0.5*(X3^3 + X3) + rnorm(n)
  Y1 <- A * (X1*(2*X1+1)^2*(X1-1)^2) + X2 + X1*X2 + 0.5*(X3^3 + X3) + rnorm(n)    # Yt
  Y <- A*Y1 + (1-A)*Y0 
  
  dat <- data.frame(X1 = X1, X2 = X2, X3 = X3, A = A, Y = Y)
  return(dat)
}


gen.dat12 <- function(n){
  
  # X = (X1, X2, X3)
  X1 <- runif(n, -0.5,0.5) # rnorm(n)   
  X2 <- sample(c(0,1,2), n, T)# rnorm(n)
  X3 <- rnorm(n)
  
  # A 
  alpha <-  c(1/2, 1/4, -1/8)*3  # rep(0.5, p)
  # alpha <-  c(1, 1, -1, -1) * 0.25
  temp <-  exp( cbind(8*X1^2, X2^2, 5*X3^2)  %*% alpha)
  # alpha <-  c(1, 1, -1) * 0.25  # rep(0.5, p)
  # temp <-  exp( cbind(X1^2, X2^2, X3^2)  %*% alpha)
  A <-  ifelse(runif(n) <=   temp / (1+temp), 1, 0)  
  
  # Y0, Y1
  Y0 <-   X2 + X1*X2 + 0.5*(X3^3 + X3) + rnorm(n)
  Y1 <- A * (cos(3*X1) * exp(X1) * log(X1+2)) + X2 + X1*X2 + 0.5*(X3^3 + X3) + rnorm(n)   
  Y <- A*Y1 + (1-A)*Y0 
  
  dat <- data.frame(X1 = X1, X2 = X2, X3 = X3, A = A, Y = Y)
  return(dat)
}



true.tau1 <- function(z){
  2*z*z
}

true.tau2 <- function(z){
  z*(2*z+1)^2*(z-1)^2
}

true.tau3 <- function(z){
  cos(3*z) * exp(z) * log(z + 2)
}

true.tau4 <- function(z){
  2*z*z
}

true.tau5 <- function(z){
  z*(2*z+1)^2*(z-1)^2
}

true.tau6 <- function(z){
  cos(3*z) * exp(z) * log(z + 2)
}

true.tau7 <- function(z){
  2*z*z
}

true.tau8 <- function(z){
  z*(2*z+1)^2*(z-1)^2
}

true.tau9 <- function(z){
  cos(3*z) * exp(z) * log(z + 2)
}

true.tau10 <- function(z){
  2*z*z
}

true.tau11 <- function(z){
  z*(2*z+1)^2*(z-1)^2
}

true.tau12 <- function(z){
  cos(3*z) * exp(z) * log(z + 2)
}

