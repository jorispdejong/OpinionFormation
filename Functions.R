# libraries
suppressPackageStartupMessages(library('pracma')) # Floating Point Relative Accuracy (eps in Matlab)

#####################################
### Compute the propensity vector ###
#####################################
propensityVector <- function(X){
  c(kappa1*(L-X[1])*exp(a1*X[1]+a2*X[2]), kappa1*(L+X[1])*exp(-a1*X[1]-a2*X[2]),
    kappa2*(L-X[2])*exp(a3*X[1]), kappa2*(L+X[2])*exp(-a3*X[1]))
}

###########################################
### Krylov algorithm for w = exp(t*A)*v ###
###########################################
# This function has been rewritten in R from the Matlab package EXPOKIT
expv <- function(t, A, v, tol, m){
  n <- nrow(A)
  anorm <- norm(A,'f') # infinity norm 
  mxrej <- 10  
  btol <- 1e-7
  gamma <- 0.9 
  delta <- 1.2
  mb <- m 
  t_out <- abs(t)
  nstep <- 0
  t_new <- 0
  t_now <- 0 
  s_error <- 0
  rndoff <- anorm*eps()
  
  k1 <- 2 
  xm <- 1/m 
  normv <- norm(v, '2')
  beta <- normv
  fact <- (((m+1)/exp(1))^(m+1))*sqrt(2*pi*(m+1))
  t_new <- (1/anorm)*((fact*tol)/(4*beta*anorm))^xm
  s <- 10^(floor(log10(t_new))-1) 
  t_new <- ceil(t_new/s)*s
  sgn <- sign(t)
  nstep <- 0
  
  w <- v
  hump <- normv
  while(t_now < t_out){
    nstep <- nstep + 1
    t_step <- min(t_out-t_now,t_new)
    V <- matrix(0,n,m+1)
    H <- matrix(0,m+2,m+2)
    
    V[,1] = (1/beta)*w
    for(j in 1:m){
      p <- A%*%V[,j]
      for(i in 1:j){
        H[i,j] <- as.numeric(t(V[,i])%*%p)
        p <- p-H[i,j]*V[,i]
      }
      s <- norm(p, '2')
      if(s < btol){
        k1 <- 0
        mb <- j
        t_step <- t_out-t_now
        break
      }
      H[j+1,j] <- s
      V[,j+1] <- as.matrix((1/s)*p)
    }
    if(k1 != 0){
      H[m+2,m+1] <- 1
      avnorm <- norm(A%*%V[,m+1], '2') 
    }
    ireject <- 0
    while(ireject <= mxrej){
      mx <- mb + k1
      F <- expm(as.matrix(sgn*t_step*H[1:mx,1:mx]))
      if(k1 == 0){
        err_loc <- btol
        break
      }
      else{
        phi1 <- abs(beta*F[m+1,1])
        phi2 <- abs(beta*F[m+2,1]) * avnorm
        if(phi1 > 10*phi2){
          err_loc <- phi2
          xm <- 1/m
        }else if(phi1 > phi2){
          err_loc <- (phi1*phi2)/(phi1-phi2)
          xm <- 1/m 
        }else{
          err_loc <- phi1
          xm <- 1/(m-1) 
        }
      }
      if(err_loc <= delta * t_step * tol){
        break 
      }
      else{
        t_step <- gamma * t_step * (t_step*tol/err_loc)^xm
        s <- 10^(floor(log10(t_step))-1)
        t_step <- ceil(t_step/s) * s
        if(ireject == mxrej){
          stop('The requested tolerance is too high')
        }
        ireject <- ireject + 1
      }
    }
    
    err_loc <- max(err_loc,rndoff)
    
    mx <- mb + max(0,k1-1)
    w <- as.matrix(V[,1:mx]) %*% as.matrix((beta*F[1:mx,1]))
    beta <- norm(w, '2')
    hump <- max(hump,beta)
    
    t_now <- t_now + t_step
    t_new <- gamma * t_step * (t_step*tol/err_loc)^xm
    s <- 10^(floor(log10(t_new))-1)
    t_new <- ceil(t_new/s) * s
    
    s_error <- s_error + err_loc
  }
  
  err <- s_error
  hump <- hump / normv
  
  return(w)
}