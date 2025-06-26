

# Match, based on X

# X: covariates used to matching. 
# A \in {0, 1}: binary treatment.
# Y: outcome. 

# matching-based estimator 
match_y1y0 <- function(X, A, Y, K = 5){
  
  n <- length(A)
  Y1 <- Y0 <- rep(NA, n)
  Y1[A==1] <- Y[A==1]
  Y0[A==0] <- Y[A==0]
  
  # matching with "manhattan" distance 
  X <- scale(X)
  distance = dist(X)
  distance = as.matrix(distance)
  
  # conduct matching
  sub0 = which(A == 0)
  sub1 = which(A == 1)
  
  for(j in sub1){
    dd = distance[j,]
    # ind1 = sub0[which.min(dd[sub0])]
    # Y0[j] <- Y0[ind1]
    ind1 = sub0[order(dd[sub0])[1:K]]
    Y0[j] <- mean(Y0[ind1]) 
  }
  
  for(j in sub0){
    dd = distance[j,]
    ind1 = sub1[order(dd[sub1])[1:K]]
    Y1[j] <- mean(Y1[ind1]) 
  }
  
  return(data.frame(Y1 = Y1, Y0 = Y0))
}


# bias-corrected matching estimator 
match_y1y0_bc <- function(X, A, Y, miu1.hat, miu0.hat, K = 5){
  
  n <- length(A)
  Y1 <- Y0 <- rep(NA, n)
  Y1[A==1] <- Y[A==1]
  Y0[A==0] <- Y[A==0]
  
  # matching with "manhattan" distance 
  X <- scale(X)
  distance = dist(X)
  distance = as.matrix(distance)
  
  # conduct matching
  sub0 = which(A == 0)
  sub1 = which(A == 1)
  
  for(j in sub1){
    dd = distance[j,]
    # ind1 = sub0[which.min(dd[sub0])]
    # Y0[j] <- Y0[ind1]
    ind1 = sub0[order(dd[sub0])[1:K]]
    Y0[j] <- miu0.hat[j] + mean(Y0[ind1] - miu0.hat[ind1]) 
  }
  
  for(j in sub0){
    dd = distance[j,]
    ind1 = sub1[order(dd[sub1])[1:K]]
    Y1[j] <- miu1.hat[j] + mean(Y1[ind1] - miu1.hat[ind1]) 
  }
  
  return(data.frame(Y1 = Y1, Y0 = Y0))
}