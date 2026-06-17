###
# Transformed linear operations
f <- function(y){
  x <- ifelse(y > 100, y, log(1+exp(y)))
  return(x)
}

finv <- function(x){
  y <- ifelse(x > 100, x, log(exp(x)-1))
}

tadd <- function(a,b){
  return(f(finv(a)+finv(b)))
}

tmult <- function(a,b){
  return(f(a*finv(b)))
}

###
# calculate tpdf
TPDF <- function(ts, thresh = 0.975, maxlag = 50){
  rt0 <- sqrt(ts^2 + ts^2)
  r00 <- sort(rt0)[floor(thresh*(length(rt0)+1))]
  tpdf <- rep(0, maxlag+1)
  
  tpdf[1] <- 2/length(which(rt0>r00)) * sum((ts^2)[which(rt0>r00)] /
                                              (rt0[which(rt0>r00)]^2))
  
  n <- length(ts)
  rt <- matrix(NA, nrow = n-1, ncol = maxlag)
  for(i in 1:maxlag){
    rt[,i] <- c(sqrt(ts^2 + stats::lag(ts, i)^2), rep(NA, i-1))
  }
  
  r0 <- rep(NA, maxlag)
  for (i in 1:maxlag){
    temp <- na.omit(rt[,i])
    if(length(temp)){ # ignore zero-length columns
      k <- floor(thresh*(length(temp)+1))
      r0[i] <- sort(temp)[k]
    } else {r0[i] <- 0}
  }
  
  for(i in 1:maxlag){
    indices <- which(rt[,i]>r0[i])
    tpdf[i+1] <- 2/length(indices) * sum((ts*stats::lag(ts, i))[indices] /
                                           (rt[indices, i]^2))
  }
  
  return(tpdf)
}

#tpdf from MA(q) thetas
MAq_tpdf <- function(theta){
  q <- length(theta)
  sigma <- rep(NA, q)
  
  theta_0 <- ts(c(1,theta))
  for(i in 0:q){
    sigma[i+1] <- sum(ifelse(theta_0>0, theta_0, 0) * 
                        stats::lag(ifelse(theta_0>0, theta_0, 0), i))/sum(theta_0^2)
  }
  return(sigma)
}


innovations <- function(tpdf, max_q =50){
  #initialize variables nu and theta_hat
  nu <- rep(NA, max_q + 1) 
  theta_hat <- matrix(0, nrow = max_q, ncol = max_q)
  
  #set nu_0 equal to tpdf(0)
  nu[1] <- tpdf[1]
  
  #compute n theta_hat's and nu_n up to n=max_q iterations
  for(n in 1:max_q){
    for(k in 0:(n-1)){
      if(k==0){
        temp <- tpdf[n+1]
      }else {
        temp <- 0
        for(j in 0:(k-1)){
          temp <- temp + theta_hat[k, k-j]*theta_hat[n, n-j]*nu[j+1]
        }
        temp <- tpdf[n-k+1] - temp
      }
      theta_hat[n, n-k] <- nu[k+1]^(-1)*temp
    }
    
    nu[n+1] <- tpdf[1]
    for(l in 0:(n-1)){
      nu[n+1] <- nu[n+1] - theta_hat[n, n-l]^2*nu[l+1]
    }
  }
  
  results <- list()
  results[[1]] <- theta_hat
  results[[2]] <- nu
  
  return(results)
}
