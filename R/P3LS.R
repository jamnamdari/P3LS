

library(gradfps)
library(fps)  # The modified version
library(RSpectra)
library(mvtnorm)
library(Matrix)
library(astsa)
library(lpSolve)
library(lattice)
library(astsa)
library(fda)
library(nhppp)
library(pracma)
library(ggplot2)
library(mnormt)

library(R.matlab)
library(spatstat)
library(compiler)
library(parallel)
library(foreach)
library(doParallel)
library(doSNOW)
library(refund)
library(fdapace)

##################################################################################

inprod_mat <- function(u1,u2,K_hat, lower, upper){
  del <- (upper-lower)/length(u1)
  value <- (del^2)*t(u1)%*%K_hat%*%u2
  return(value[1,1])
}



GS <- function(u_list,K_hat, lower, upper){

  u1_1 <- u_list[[1]]
  u1_1_norm <- sqrt(inprod_mat(u1_1,u1_1,K_hat, lower, upper))
  u1 <- u1_1/u1_1_norm

  u2_1 <- u_list[[2]]
  u2_2 <- u2_1 - inprod_mat(u2_1,u1,K_hat, lower, upper)*u1
  u2_2_norm <- sqrt(inprod_mat(u2_2,u2_2,K_hat, lower, upper))
  u2 <-  u2_2/u2_2_norm

  u3_1 <- u_list[[3]]
  u3_2 <- u3_1 - inprod_mat(u3_1,u1,K_hat, lower, upper)*u1
  u3_3 <- u3_2 - inprod_mat(u3_2,u2,K_hat, lower, upper)*u2
  u3_3_norm <- sqrt(inprod_mat(u3_3,u3_3,K_hat, lower, upper))
  u3 <- u3_3/u3_3_norm

  u4_1 <- u_list[[4]]
  u4_2 <- u4_1 - inprod_mat(u4_1,u1,K_hat, lower, upper)*u1
  u4_3 <- u4_2 - inprod_mat(u4_2,u2,K_hat, lower, upper)*u2
  u4_4 <- u4_3 - inprod_mat(u4_3,u3,K_hat, lower, upper)*u3
  u4_4_norm <- sqrt(inprod_mat(u4_4,u4_4,K_hat, lower, upper))
  u4 <- u4_4/u4_4_norm

  u5_1 <- u_list[[5]]
  u5_2 <- u5_1 - inprod_mat(u5_1,u1,K_hat, lower, upper)*u1
  u5_3 <- u5_2 - inprod_mat(u5_2,u2,K_hat, lower, upper)*u2
  u5_4 <- u5_3 - inprod_mat(u5_3,u3,K_hat, lower, upper)*u3
  u5_5 <- u5_4 - inprod_mat(u5_4,u4,K_hat, lower, upper)*u4
  u5_5_norm <- sqrt(inprod_mat(u5_5,u5_5,K_hat, lower, upper))
  u5 <- u5_5/u5_5_norm

  u6_1 <- u_list[[6]]
  u6_2 <- u6_1 - inprod_mat(u6_1,u1,K_hat, lower, upper)*u1
  u6_3 <- u6_2 - inprod_mat(u6_2,u2,K_hat, lower, upper)*u2
  u6_4 <- u6_3 - inprod_mat(u6_3,u3,K_hat, lower, upper)*u3
  u6_5 <- u6_4 - inprod_mat(u6_4,u4,K_hat, lower, upper)*u4
  u6_6 <- u6_5 - inprod_mat(u6_5,u5,K_hat, lower, upper)*u5
  u6_6_norm <- sqrt(inprod_mat(u6_6,u6_6,K_hat, lower, upper))
  u6 <- u6_6/u6_6_norm

  u7_1 <- u_list[[7]]
  u7_2 <- u7_1 - inprod_mat(u7_1,u1,K_hat, lower, upper)*u1
  u7_3 <- u7_2 - inprod_mat(u7_2,u2,K_hat, lower, upper)*u2
  u7_4 <- u7_3 - inprod_mat(u7_3,u3,K_hat, lower, upper)*u3
  u7_5 <- u7_4 - inprod_mat(u7_4,u4,K_hat, lower, upper)*u4
  u7_6 <- u7_5 - inprod_mat(u7_5,u5,K_hat, lower, upper)*u5
  u7_7 <- u7_6 - inprod_mat(u7_6,u6,K_hat, lower, upper)*u6
  u7_7_norm <- sqrt(inprod_mat(u7_7,u7_7,K_hat, lower, upper))
  u7 <- u7_7/u7_7_norm

  u8_1 <- u_list[[8]]
  u8_2 <- u8_1 - inprod_mat(u8_1,u1,K_hat, lower, upper)*u1
  u8_3 <- u8_2 - inprod_mat(u8_2,u2,K_hat, lower, upper)*u2
  u8_4 <- u8_3 - inprod_mat(u8_3,u3,K_hat, lower, upper)*u3
  u8_5 <- u8_4 - inprod_mat(u8_4,u4,K_hat, lower, upper)*u4
  u8_6 <- u8_5 - inprod_mat(u8_5,u5,K_hat, lower, upper)*u5
  u8_7 <- u8_6 - inprod_mat(u8_6,u6,K_hat, lower, upper)*u6
  u8_8 <- u8_7 - inprod_mat(u8_7,u7,K_hat, lower, upper)*u7
  u8_8_norm <- sqrt(inprod_mat(u8_8,u8_8,K_hat, lower, upper))
  u8 <- u8_8/u8_8_norm

  u9_1 <- u_list[[9]]
  u9_2 <- u9_1 - inprod_mat(u9_1,u1,K_hat, lower, upper)*u1
  u9_3 <- u9_2 - inprod_mat(u9_2,u2,K_hat, lower, upper)*u2
  u9_4 <- u9_3 - inprod_mat(u9_3,u3,K_hat, lower, upper)*u3
  u9_5 <- u9_4 - inprod_mat(u9_4,u4,K_hat, lower, upper)*u4
  u9_6 <- u9_5 - inprod_mat(u9_5,u5,K_hat, lower, upper)*u5
  u9_7 <- u9_6 - inprod_mat(u9_6,u6,K_hat, lower, upper)*u6
  u9_8 <- u9_7 - inprod_mat(u9_7,u7,K_hat, lower, upper)*u7
  u9_9 <- u9_8 - inprod_mat(u9_8,u8,K_hat, lower, upper)*u8
  u9_9_norm <- sqrt(inprod_mat(u9_9,u9_9,K_hat, lower, upper))
  u9 <- u9_9/u9_9_norm

  u10_1 <- u_list[[10]]
  u10_2 <- u10_1 - inprod_mat(u10_1,u1,K_hat, lower, upper)*u1
  u10_3 <- u10_2 - inprod_mat(u10_2,u2,K_hat, lower, upper)*u2
  u10_4 <- u10_3 - inprod_mat(u10_3,u3,K_hat, lower, upper)*u3
  u10_5 <- u10_4 - inprod_mat(u10_4,u4,K_hat, lower, upper)*u4
  u10_6 <- u10_5 - inprod_mat(u10_5,u5,K_hat, lower, upper)*u5
  u10_7 <- u10_6 - inprod_mat(u10_6,u6,K_hat, lower, upper)*u6
  u10_8 <- u10_7 - inprod_mat(u10_7,u7,K_hat, lower, upper)*u7
  u10_9 <- u10_8 - inprod_mat(u10_8,u8,K_hat, lower, upper)*u8
  u10_10 <- u10_9 - inprod_mat(u10_9,u9,K_hat, lower, upper)*u9
  u10_10_norm <- sqrt(inprod_mat(u10_10,u10_10,K_hat, lower, upper))
  u10 <- u10_10/u10_10_norm

  return(list(u1,u2,u3,u4,u5,u6,u7,u8,u9,u10))
}


### **Covariance estimation by point process method**

Cov_estimator <- function(PROCESS, lbd, ubd, bwd, ngrid, kern="epanechnikov"){
  lower <- lbd
  upper <- ubd

  n <- length(PROCESS)
  grids <- seq(lbd,ubd,l=ngrid)
  edge <- pkernel((grids-lbd)/bwd,kernel = kern)-pkernel((grids-ubd)/bwd,kernel = kern,sd=sqrt(1/5))
  edge <- outer(edge,edge,FUN="*")

  Kh <- function(t){
    dkernel(t/bwd,kernel = kern,sd=sqrt(1/5))/bwd
  }
  A2 <- A <- B <- C <- D <- matrix(0,ngrid,ngrid)

  #start_time <- Sys.time()
  m <- 1
  A2.out <- list()
  A2.out.y <- vector(mode = "list", length = m)
  for(i in 1:n){
    tmp.A2 <- matrix(0,ngrid,ngrid)
    for(j in 1:m){
      process <- PROCESS[[i]]
      if(!is.null(process)){
        cp <- crosspairs(ppp(x=process,y=process*0, window = owin(c(lbd,ubd),c(lbd,ubd))),ppp(x=grids,y=grids*0, window = owin(c(lbd,ubd),c(lbd,ubd))),rmax = bwd,what = "ijd")
        tmp <- Matrix::sparseMatrix(i = cp$i, j = cp$j, x = Kh(cp$d), dims = c(length(process), ngrid))
        tmp1 <- apply(tmp,2,sum)
        tmp2 <- outer(tmp1,tmp1,FUN="*")
        A <- A+(tmp2-Matrix::t(tmp)%*%tmp)
        A2 <- A2 + tmp2
        tmp.A2 <- tmp.A2 + tmp2
        if(is.null(A2.out.y[[j]])) A2.out.y[[j]] <- tmp2 else
          A2.out.y[[j]] <- A2.out.y[[j]] +tmp2
      }
    }
    A2.out[[i]] <- tmp.A2
    if(i%%50==0) print(i)
  }


  nlarge <- 10^6
  process.y <- rho.out.y <- C.out <- list()
  process_l <- PROCESS
  for(j in 1:m){
    process.y[[j]] <- process.tmp <- unlist(process_l)
    if(length(process.tmp)>nlarge){
      grid.large <- seq(lbd,ubd,l=nlarge+1)
      grid.mid <- (grid.large[-1]+grid.large[-(nlarge+1)])/2
      process <- grid.mid
      freq <- c(table(cut(process.tmp,breaks=grid.large)))
    }else{
      process <- process.tmp
      freq <- rep(1,length(process.tmp))
    }
    cp <- crosspairs(ppp(x=process,y=process*0, window = owin(c(lbd,ubd),c(lbd,ubd))),ppp(x=grids,y=grids*0, window = owin(c(lbd,ubd),c(lbd,ubd))),rmax = bwd,what = "ijd")
    tmp <- Matrix::sparseMatrix(i = cp$i, j = cp$j, x = Kh(cp$d)*freq[cp$i], dims = c(length(process), ngrid))
    tmp1 <- apply(tmp,2,sum)
    tmp2 <- outer(tmp1,tmp1,FUN="*")
    C <- C + tmp2
    C.out[[j]] <- (tmp2-A2.out.y[[j]])/(n*(n-1)*edge)
    rho.out.y[[j]] <- c(tmp1)/sqrt(diag(edge))/n
  }


  C <- C-A2

  A <- as.matrix(A/(n*m*edge))
  C <- as.matrix(C/((n*m*(n-1))*edge))

  R_X <- log(A/C)

  #end_time <- Sys.time()
  #time[l] <- difftime(end_time,start_time,units="secs")
  return(R_X)

}

########################

#########################################
## Using Kernel Method
#########################################
###############################
## Set up
###############################


Ker <- function(s,u,h){
  v <- (s-u)/h
  ifelse(abs(v)<1,(1/h)*(3/4)*(1-(v)^2),0)
}
edge_c <- function(s,h,lower,upper) {
  integrate(function(x) Ker(s,x,h), lower = lower, upper = upper)$value
}

#' @export
lambda_fun <- function(PP, s,h){
  sum(Ker(s,PP,h))/edge_c(s,h, lower, upper)
}


#' @export
PLS_Kernel <- function(PPP_obs, PPP_test = NULL, y, y_test=NULL, h, T, lbd, ubd){
  lower <- lbd
  upper <- ubd
  Ker <- function(s,u,h){
    v <- (s-u)/h
    ifelse(abs(v)<1,(1/h)*(3/4)*(1-(v)^2),0)
  }
  edge_c <- function(s,h,lower=lbd,upper=ubd) {
    integrate(function(x) Ker(s,x,h), lower = lower, upper = upper)$value
  }

  lambda_fun <- function(PP, s,h){
    sum(Ker(s,PP,h))/edge_c(s,h, lower, upper)
  }



  n_obs <- length(PPP_obs)
  if(!is.null(PPP_test)){
    n_test <- length(PPP_test)
  }
  n <- n_obs
  dd <- 10
  b_hat_kr  <- matrix(0,nrow = T, ncol=dd)
  y_hat_kr  <- matrix(0,nrow = T, ncol=dd)
  y_kr_err_l2  <- rep(0,dd)
  yhat_test_K  <- matrix(0,nrow = T, ncol=dd)
  y_test_K_err_l2  <- rep(0,dd)

  yc <- y - mean(y)
  y_bar <- mean(y)


  del <- (upper-lower)/T
  T_seq <- seq(lbd,ubd,length.out=T)

  X_k <- matrix(unlist(lapply(PPP_obs, function(PP) sapply(T_seq, function(s) log(lambda_fun(PP,s,h))))), nrow = n_obs, byrow = TRUE)
  Xkc <- t(t(X_k) - colMeans(X_k))

  Kk_hat <- (1/n_obs)*t(Xkc)%*%Xkc
  K_b <- (1/n_obs)*t(Xkc)%*%yc
  K2_b <- (del)*Kk_hat%*%K_b
  K3_b <- ((del))*Kk_hat%*%K2_b
  K4_b <- ((del))*Kk_hat%*%K3_b
  K5_b <- ((del))*Kk_hat%*%K4_b
  K6_b <- ((del))*Kk_hat%*%K5_b
  K7_b <- ((del))*Kk_hat%*%K6_b
  K8_b <- ((del))*Kk_hat%*%K7_b
  K9_b <- ((del))*Kk_hat%*%K8_b
  K10_b <- ((del))*Kk_hat%*%K9_b

  U <- list(K_b,K2_b,K3_b,K4_b,K5_b,K6_b,K7_b,K8_b,K9_b,K10_b)
  V <- GS(U,Kk_hat, lower, upper)

  AIC_K <- rep(0,dd)
  BIC_K <- rep(0,dd)
  SSE_K <- rep(0,dd)

  if(!is.null(PPP_test)){
    X_k_test <- matrix(unlist(lapply(PPP_test, function(PP) sapply(T_seq, function(s) log(lambda_fun(PP,s,h))))), nrow = n_test, byrow = TRUE)
    Xkc_test <- t(t(X_k_test) - colMeans(X_k_test))
  }



  p=1
  orth_b <- cbind(V[[1]])
  del <- (upper-lower)/T
  A <- del*Xkc%*%orth_b
  w <- inv(t(A)%*%A)%*%t(A)%*%yc
  b_hat_kr[,p] <- w[1]*V[[1]]
  y_hat_kr[,p] <- w[1]*A[,1]
  SSE_K[p] <- sum((yc-y_hat_kr[,p])^2)/n
  AIC_K[p] <- n*log(SSE_K[p]) + 2*p
  BIC_K[p] <- n*log(SSE_K[p]) + p*log(n)
  y_kr_err_l2[p] <- sqrt(SSE_K[p])

  if(!is.null(PPP_test)){
    yhat_test_K[,p] <- y_bar + del*(t(Xkc_test)%*%b_hat_kr[,p])
    y_test_K_err_l2[p] <- sqrt(sum((y_test-yhat_test_K[,p])^2)/n)
  }

  p=2
  orth_b <- cbind(V[[1]],V[[2]])
  del <- (upper-lower)/T
  A <- del*Xkc%*%orth_b
  w <- inv(t(A)%*%A)%*%t(A)%*%yc
  b_hat_kr[,p] <- w[1]*V[[1]]+w[2]*V[[2]]
  y_hat_kr[,p] <- w[1]*A[,1]+w[2]*A[,2]
  SSE_K[p] <- sum((yc-y_hat_kr[,p])^2)/n
  AIC_K[p] <- n*log(SSE_K[p]) + 2*p
  BIC_K[p] <- n*log(SSE_K[p]) + p*log(n)
  y_kr_err_l2[p] <- sqrt(SSE_K[p])

  if(!is.null(PPP_test)){
    yhat_test_K[,p] <- y_bar + del*(t(Xkc_test)%*%b_hat_kr[,p])
    y_test_K_err_l2[p] <- sqrt(sum((y_test-yhat_test_K[,p])^2)/n)
  }

  p=3
  orth_b <- cbind(V[[1]],V[[2]],V[[3]])
  del <- (upper-lower)/T
  A <- del*Xkc%*%orth_b
  w <- inv(t(A)%*%A)%*%t(A)%*%yc
  b_hat_kr[,p] <- w[1]*V[[1]]+w[2]*V[[2]]+w[3]*V[[3]]
  y_hat_kr[,p] <- w[1]*A[,1]+w[2]*A[,2]+w[3]*A[,3]
  SSE_K[p] <- sum((yc-y_hat_kr[,p])^2)/n
  AIC_K[p] <- n*log(SSE_K[p]) + 2*p
  BIC_K[p] <- n*log(SSE_K[p]) + p*log(n)
  y_kr_err_l2[p] <- sqrt(SSE_K[p])

  if(!is.null(PPP_test)){
    yhat_test_K[,p] <- y_bar + del*(t(Xkc_test)%*%b_hat_kr[,p])
    y_test_K_err_l2[p] <- sqrt(sum((y_test-yhat_test_K[,p])^2)/n)
  }

  p=4
  orth_b <- cbind(V[[1]],V[[2]],V[[3]],V[[4]])
  del <- (upper-lower)/T
  A <- del*Xkc%*%orth_b
  w <- inv(t(A)%*%A)%*%t(A)%*%yc
  b_hat_kr[,p] <- w[1]*V[[1]]+w[2]*V[[2]]+w[3]*V[[3]]+w[4]*V[[4]]
  y_hat_kr[,p] <- w[1]*A[,1]+w[2]*A[,2]+w[3]*A[,3]+w[4]*A[,4]
  SSE_K[p] <- sum((yc-y_hat_kr[,p])^2)/n
  AIC_K[p] <- n*log(SSE_K[p]) + 2*p
  BIC_K[p] <- n*log(SSE_K[p]) + p*log(n)
  y_kr_err_l2[p] <- sqrt(SSE_K[p])

  if(!is.null(PPP_test)){
    yhat_test_K[,p] <- y_bar + del*(t(Xkc_test)%*%b_hat_kr[,p])
    y_test_K_err_l2[p] <- sqrt(sum((y_test-yhat_test_K[,p])^2)/n)
  }

  p=5
  orth_b <- cbind(V[[1]],V[[2]],V[[3]],V[[4]],V[[5]])
  del <- (upper-lower)/T
  A <- del*Xkc%*%orth_b
  w <- ginv(t(A)%*%A)%*%t(A)%*%yc
  b_hat_kr[,p] <- w[1]*V[[1]]+w[2]*V[[2]]+w[3]*V[[3]]+w[4]*V[[4]]+w[5]*V[[5]]
  y_hat_kr[,p] <- w[1]*A[,1]+w[2]*A[,2]+w[3]*A[,3]+w[4]*A[,4]+w[5]*A[,5]
  SSE_K[p] <- sum((yc-y_hat_kr[,p])^2)/n
  AIC_K[p] <- n*log(SSE_K[p]) + 2*p
  BIC_K[p] <- n*log(SSE_K[p]) + p*log(n)
  y_kr_err_l2[p] <- sqrt(SSE_K[p])

  if(!is.null(PPP_test)){
    yhat_test_K[,p] <- y_bar + del*(t(Xkc_test)%*%b_hat_kr[,p])
    y_test_K_err_l2[p] <- sqrt(sum((y_test-yhat_test_K[,p])^2)/n)
  }

  p=6
  orth_b <- cbind(V[[1]],V[[2]],V[[3]],V[[4]],V[[5]],V[[6]])
  del <- (upper-lower)/T
  A <- del*Xkc%*%orth_b
  w <- ginv(t(A)%*%A)%*%t(A)%*%yc
  b_hat_kr[,p] <- w[1]*V[[1]]+w[2]*V[[2]]+w[3]*V[[3]]+w[4]*V[[4]]+w[5]*V[[5]]+w[6]*V[[6]]
  y_hat_kr[,p] <- w[1]*A[,1]+w[2]*A[,2]+w[3]*A[,3]+w[4]*A[,4]+w[5]*A[,5]+w[6]*A[,6]
  SSE_K[p] <- sum((yc-y_hat_kr[,p])^2)/n
  AIC_K[p] <- n*log(SSE_K[p]) + 2*p
  BIC_K[p] <- n*log(SSE_K[p]) + p*log(n)
  y_kr_err_l2[p] <- sqrt(SSE_K[p])

  if(!is.null(PPP_test)){
    yhat_test_K[,p] <- y_bar + del*(t(Xkc_test)%*%b_hat_kr[,p])
    y_test_K_err_l2[p] <- sqrt(sum((y_test-yhat_test_K[,p])^2)/n)
  }

  p=7
  orth_b <- cbind(V[[1]],V[[2]],V[[3]],V[[4]],V[[5]],V[[6]],V[[7]])
  del <- (upper-lower)/T
  A <- del*Xkc%*%orth_b
  w <- ginv(t(A)%*%A)%*%t(A)%*%yc
  b_hat_kr[,p] <- w[1]*V[[1]]+w[2]*V[[2]]+w[3]*V[[3]]+w[4]*V[[4]]+w[5]*V[[5]]+w[6]*V[[6]]+w[7]*V[[7]]
  y_hat_kr[,p] <- w[1]*A[,1]+w[2]*A[,2]+w[3]*A[,3]+w[4]*A[,4]+w[5]*A[,5]+w[6]*A[,6]+w[7]*A[,7]
  SSE_K[p] <- sum((yc-y_hat_kr[,p])^2)/n
  AIC_K[p] <- n*log(SSE_K[p]) + 2*p
  BIC_K[p] <- n*log(SSE_K[p]) + p*log(n)
  y_kr_err_l2[p] <- sqrt(SSE_K[p])

  if(!is.null(PPP_test)){
    yhat_test_K[,p] <- y_bar + del*(t(Xkc_test)%*%b_hat_kr[,p])
    y_test_K_err_l2[p] <- sqrt(sum((y_test-yhat_test_K[,p])^2)/n)
  }

  p=8
  orth_b <- cbind(V[[1]],V[[2]],V[[3]],V[[4]],V[[5]],V[[6]],V[[7]],V[[8]])
  del <- (upper-lower)/T
  A <- del*Xkc%*%orth_b
  w <- ginv(t(A)%*%A)%*%t(A)%*%yc
  b_hat_kr[,p] <- w[1]*V[[1]]+w[2]*V[[2]]+w[3]*V[[3]]+w[4]*V[[4]]+w[5]*V[[5]]+w[6]*V[[6]]+w[7]*V[[7]]+w[8]*V[[8]]
  y_hat_kr[,p] <- w[1]*A[,1]+w[2]*A[,2]+w[3]*A[,3]+w[4]*A[,4]+w[5]*A[,5]+w[6]*A[,6]+w[7]*A[,7]+w[8]*A[,8]
  SSE_K[p] <- sum((yc-y_hat_kr[,p])^2)/n
  AIC_K[p] <- n*log(SSE_K[p]) + 2*p
  BIC_K[p] <- n*log(SSE_K[p]) + p*log(n)
  y_kr_err_l2[p] <- sqrt(SSE_K[p])

  if(!is.null(PPP_test)){
    yhat_test_K[,p] <- y_bar + del*(t(Xkc_test)%*%b_hat_kr[,p])
    y_test_K_err_l2[p] <- sqrt(sum((y_test-yhat_test_K[,p])^2)/n)
  }

  p=9
  orth_b <- cbind(V[[1]],V[[2]],V[[3]],V[[4]],V[[5]],V[[6]],V[[7]],V[[8]],V[[9]])
  del <- (upper-lower)/T
  A <- del*Xkc%*%orth_b
  w <- ginv(t(A)%*%A)%*%t(A)%*%yc
  b_hat_kr[,p] <- w[1]*V[[1]]+w[2]*V[[2]]+w[3]*V[[3]]+w[4]*V[[4]]+w[5]*V[[5]]+w[6]*V[[6]]+w[7]*V[[7]]+w[8]*V[[8]]+w[9]*V[[9]]
  y_hat_kr[,p] <- w[1]*A[,1]+w[2]*A[,2]+w[3]*A[,3]+w[4]*A[,4]+w[5]*A[,5]+w[6]*A[,6]+w[7]*A[,7]+w[8]*A[,8]+w[9]*A[,9]
  SSE_K[p] <- sum((yc-y_hat_kr[,p])^2)/n
  AIC_K[p] <- n*log(SSE_K[p]) + 2*p
  BIC_K[p] <- n*log(SSE_K[p]) + p*log(n)
  y_kr_err_l2[p] <- sqrt(SSE_K[p])

  if(!is.null(PPP_test)){
    yhat_test_K[,p] <- y_bar + del*(t(Xkc_test)%*%b_hat_kr[,p])
    y_test_K_err_l2[p] <- sqrt(sum((y_test-yhat_test_K[,p])^2)/n)
  }

  p=10
  orth_b <- cbind(V[[1]],V[[2]],V[[3]],V[[4]],V[[5]],V[[6]],V[[7]],V[[8]],V[[9]],V[[10]])
  del <- (upper-lower)/T
  A <- del*Xkc%*%orth_b
  w <- ginv(t(A)%*%A)%*%t(A)%*%yc
  b_hat_kr[,p] <- w[1]*V[[1]]+w[2]*V[[2]]+w[3]*V[[3]]+w[4]*V[[4]]+w[5]*V[[5]]+w[6]*V[[6]]+w[7]*V[[7]]+w[8]*V[[8]]+w[9]*V[[9]]+w[10]*V[[10]]
  y_hat_kr[,p] <- w[1]*A[,1]+w[2]*A[,2]+w[3]*A[,3]+w[4]*A[,4]+w[5]*A[,5]+w[6]*A[,6]+w[7]*A[,7]+w[8]*A[,8]+w[9]*A[,9]+w[10]*A[,10]
  SSE_K[p] <- sum((yc-y_hat_kr[,p])^2)/n
  AIC_K[p] <- n*log(SSE_K[p]) + 2*p
  BIC_K[p] <- n*log(SSE_K[p]) + p*log(n)
  y_kr_err_l2[p] <- sqrt(SSE_K[p])

  if(!is.null(PPP_test)){
    yhat_test_K[,p] <- y_bar + del*(t(Xkc_test)%*%b_hat_kr[,p])
    y_test_K_err_l2[p] <- sqrt(sum((y_test-yhat_test_K[,p])^2)/n)
  }

  nbasis_K <- order(BIC_K)[1]

  if(is.null(PPP_test)){
    yhat_test_K <- NULL
    y_test_K_err_l2 <- NULL
  }
  return(list(b_hat = b_hat_kr, y_hat = y_hat_kr, BIC = BIC_K, y_hat_test = yhat_test_K, PLS_basis = orth_b, rMSPE_test = y_test_K_err_l2))

}


#########################################
## using GLM
#########################################

#########################################
## using GLM + PPP new
#########################################
###############################
## Set up
###############################

#' @export
P3LS <- function(PPP_obs, PPP_test = NULL, y, y_test=NULL, h, p=10, q, T, lbd, ubd, nb = 100){
  p_pls <- p
  lower <- lbd
  upper <- ubd

  n_obs <- length(PPP_obs)
  if(!is.null(PPP_test)){
    n_test <- length(PPP_test)
  }
  n <- n_obs

  dd <- 10
  b_hat_pp  <- matrix(0,nrow = T, ncol=dd)
  y_hat_pp  <- matrix(0,nrow = n_obs, ncol=dd)
  y_pp_err_l2  <- rep(0,dd)
  yhat_test_PP  <- matrix(0,nrow = n_test, ncol=dd)
  y_test_PP_err_l2  <- rep(0,dd)


  AIC_PP <- rep(0,dd)
  BIC_PP <- rep(0,dd)
  SSE_PP <- rep(0,dd)


  del <- (ubd-lbd)/T
  T_seq <- seq(lbd,ubd,length.out=T)

  yc <- y - mean(y)
  y_bar <- mean(y)

  K_X1 <- Cov_estimator(PPP_obs, lbd=lbd, ubd=ubd, bwd=h, ngrid=100)
  K_hat_pp1_e <- eigen(K_X1)
  neval_inds <- K_hat_pp1_e$values < 0
  pevals <- K_hat_pp1_e$values
  pevals[neval_inds] <- 0
  K_X <- K_hat_pp1_e$vectors%*%diag(pevals)%*%t(K_hat_pp1_e$vectors)

  #q <- 3+2+5+5 # number of eigenfunctions
  ef <- t(K_hat_pp1_e$vectors[,1:q])/sqrt(del) # eigen functions evaluated at grid points
  #nb <- 100  ## number of bins

  ## Bin construction
  ends <- seq(lbd,ubd,length.out=nb)
  bins0 <- matrix(c(0,ends,ends,0), nrow=nb+1, ncol = 2, byrow = FALSE)
  bins <- bins0[-c(1,nb+1),]
  ## Bin length
  b_l <- bins[1,2]-bins[1,1]

  # phis used in the GLM
  inds <-  t(apply(bins, MARGIN = 1, FUN = function(t) (T_seq>t[1] & T_seq <= t[2]) ) )
  temp <- 1:T
  PHIsm <-  t(apply(inds, MARGIN = 1, FUN = function(t) ef[,floor(median(temp[t]))] ) )

  #bin_counts <- matrix(0, nrow = n_obs, ncol=nb)
  X_k <- matrix(0,nrow = n_obs, ncol = T)
  for(ell in 1:n_obs){
    ## Point process data
    PP_ts <- PPP_obs[[ell]]

    ## Number of data points (events) in each bin
    bin_counts <- apply(bins, MARGIN = 1, FUN = function(t) sum(PP_ts>t[1] & PP_ts <= t[2]))/b_l
    Dm <- data.frame(counts = bin_counts, phi = PHIsm)
    Pm_model <- glm(counts~.+1  , family = poisson(link = "log"), data = Dm)
    X_k[ell,] <- colSums(diag(Pm_model$coefficients)%*%rbind(rep(1,nb),ef))
  }
  Xkc <- t(t(X_k) - colMeans(X_k))

  # test set
  if(!is.null(PPP_test)){
    X_k_test <- matrix(0,nrow = n_test, ncol = T)
    for(ell in 1:n_test){
      ## Point process data
      PP_ts <- PPP_test[[ell]]

      ## Number of data points (events) in each bin
      bin_counts <- apply(bins, MARGIN = 1, FUN = function(t) sum(PP_ts>t[1] & PP_ts <= t[2]))/b_l
      Dm <- data.frame(counts = bin_counts, phi = PHIsm)
      Pm_model <- glm(counts~.+1  , family = poisson(link = "log"), data = Dm)
      X_k_test[ell,] <- colSums(diag(Pm_model$coefficients)%*%rbind(rep(1,nb),ef))
    }
    X_k_test_c <- t(t(X_k_test) - colMeans(X_k_test))
  }

  K_hat_e <- K_hat_pp1_e
  K_b <- (1/n_obs)*t(Xkc)%*%yc
  K2_b <- (del)*K_X%*%K_b
  K3_b <- ((del))*K_X%*%K2_b
  K4_b <- ((del))*K_X%*%K3_b
  K5_b <- ((del))*K_X%*%K4_b
  K6_b <- ((del))*K_X%*%K5_b
  K7_b <- ((del))*K_X%*%K6_b
  K8_b <- ((del))*K_X%*%K7_b
  K9_b <- ((del))*K_X%*%K8_b
  K10_b <- ((del))*K_X%*%K9_b

  U <- list(K_b,K2_b,K3_b,K4_b,K5_b,K6_b,K7_b,K8_b,K9_b,K10_b)
  V <- GS(U,K_X, lower, upper)


  p=1
  orth_b <- cbind(V[[1]])
  del <- (upper-lower)/T
  A <- del*Xkc%*%orth_b
  w <- inv(t(A)%*%A)%*%t(A)%*%yc
  b_hat_pp[,p] <- w[1]*V[[1]]
  y_hat_pp[,p] <- w[1]*A[,1]
  SSE_PP[p] <- sum((yc-y_hat_pp[,p])^2)/n
  AIC_PP[p] <- n*log(SSE_PP[p]) + 2*p
  BIC_PP[p] <- n*log(SSE_PP[p]) + p*log(n)
  y_pp_err_l2[p] <- sqrt(SSE_PP[p])

  if(!is.null(PPP_test)){
    yhat_test_PP[,p] <- y_bar + del*(X_k_test_c%*%b_hat_pp[,p])
    y_test_PP_err_l2[p] <- sqrt(sum((y_test-yhat_test_PP[,p])^2)/n)
  }

  p=2
  orth_b <- cbind(V[[1]],V[[2]])
  del <- (upper-lower)/T
  A <- del*Xkc%*%orth_b
  w <- inv(t(A)%*%A)%*%t(A)%*%yc
  b_hat_pp[,p] <- w[1]*V[[1]]+w[2]*V[[2]]
  y_hat_pp[,p] <- w[1]*A[,1]+w[2]*A[,2]
  SSE_PP[p] <- sum((yc-y_hat_pp[,p])^2)/n
  AIC_PP[p] <- n*log(SSE_PP[p]) + 2*p
  BIC_PP[p] <- n*log(SSE_PP[p]) + p*log(n)
  y_pp_err_l2[p] <- sqrt(SSE_PP[p])

  if(!is.null(PPP_test)){
    yhat_test_PP[,p] <- y_bar + del*(X_k_test_c%*%b_hat_pp[,p])
    y_test_PP_err_l2[p] <- sqrt(sum((y_test-yhat_test_PP[,p])^2)/n)
  }

  if(p_pls <= 2){
    nbasis_PP <- order(BIC_PP[1:p])[1]
    return(list(b_hat = b_hat_pp[,1:p], y_hat = y_hat_pp[,1:p], BIC = BIC_PP[1:p], y_hat_test = yhat_test_PP[,1:p], basis = orth_b, rMSPE_test = y_test_PP_err_l2[1:p]))

  }

  p=3
  orth_b <- cbind(V[[1]],V[[2]],V[[3]])
  del <- (upper-lower)/T
  A <- del*Xkc%*%orth_b
  w <- inv(t(A)%*%A)%*%t(A)%*%yc
  b_hat_pp[,p] <- w[1]*V[[1]]+w[2]*V[[2]]+w[3]*V[[3]]
  y_hat_pp[,p] <- w[1]*A[,1]+w[2]*A[,2]+w[3]*A[,3]
  SSE_PP[p] <- sum((yc-y_hat_pp[,p])^2)/n
  AIC_PP[p] <- n*log(SSE_PP[p]) + 2*p
  BIC_PP[p] <- n*log(SSE_PP[p]) + p*log(n)
  y_pp_err_l2[p] <- sqrt(SSE_PP[p])

  if(!is.null(PPP_test)){
    yhat_test_PP[,p] <- y_bar + del*(X_k_test_c%*%b_hat_pp[,p])
    y_test_PP_err_l2[p] <- sqrt(sum((y_test-yhat_test_PP[,p])^2)/n)
  }

  if(p_pls == p){
    nbasis_PP <- order(BIC_PP[1:p])[1]
    return(list(b_hat = b_hat_pp[,1:p], y_hat = y_hat_pp[,1:p], BIC = BIC_PP[1:p], y_hat_test = yhat_test_PP[,1:p], basis = orth_b, rMSPE_test = y_test_PP_err_l2[1:p]))

  }

  p=4
  orth_b <- cbind(V[[1]],V[[2]],V[[3]],V[[4]])
  del <- (upper-lower)/T
  A <- del*Xkc%*%orth_b
  w <- inv(t(A)%*%A)%*%t(A)%*%yc
  b_hat_pp[,p] <- w[1]*V[[1]]+w[2]*V[[2]]+w[3]*V[[3]]+w[4]*V[[4]]
  y_hat_pp[,p] <- w[1]*A[,1]+w[2]*A[,2]+w[3]*A[,3]+w[4]*A[,4]
  SSE_PP[p] <- sum((yc-y_hat_pp[,p])^2)/n
  AIC_PP[p] <- n*log(SSE_PP[p]) + 2*p
  BIC_PP[p] <- n*log(SSE_PP[p]) + p*log(n)
  y_pp_err_l2[p] <- sqrt(SSE_PP[p])

  if(!is.null(PPP_test)){
    yhat_test_PP[,p] <- y_bar + del*(X_k_test_c%*%b_hat_pp[,p])
    y_test_PP_err_l2[p] <- sqrt(sum((y_test-yhat_test_PP[,p])^2)/n)
  }

  if(p_pls == p){
    nbasis_PP <- order(BIC_PP[1:p])[1]
    return(list(b_hat = b_hat_pp[,1:p], y_hat = y_hat_pp[,1:p], BIC = BIC_PP[1:p], y_hat_test = yhat_test_PP[,1:p], basis = orth_b, rMSPE_test = y_test_PP_err_l2[1:p]))

  }

  p=5
  orth_b <- cbind(V[[1]],V[[2]],V[[3]],V[[4]],V[[5]])
  del <- (upper-lower)/T
  A <- del*Xkc%*%orth_b
  w <- ginv(t(A)%*%A)%*%t(A)%*%yc
  b_hat_pp[,p] <- w[1]*V[[1]]+w[2]*V[[2]]+w[3]*V[[3]]+w[4]*V[[4]]+w[5]*V[[5]]
  y_hat_pp[,p] <- w[1]*A[,1]+w[2]*A[,2]+w[3]*A[,3]+w[4]*A[,4]+w[5]*A[,5]
  SSE_PP[p] <- sum((yc-y_hat_pp[,p])^2)/n
  AIC_PP[p] <- n*log(SSE_PP[p]) + 2*p
  BIC_PP[p] <- n*log(SSE_PP[p]) + p*log(n)
  y_pp_err_l2[p] <- sqrt(SSE_PP[p])

  if(!is.null(PPP_test)){
    yhat_test_PP[,p] <- y_bar + del*(X_k_test_c%*%b_hat_pp[,p])
    y_test_PP_err_l2[p] <- sqrt(sum((y_test-yhat_test_PP[,p])^2)/n)
  }

  if(p_pls == p){
    nbasis_PP <- order(BIC_PP[1:p])[1]
    return(list(b_hat = b_hat_pp[,1:p], y_hat = y_hat_pp[,1:p], BIC = BIC_PP[1:p], y_hat_test = yhat_test_PP[,1:p], basis = orth_b, rMSPE_test = y_test_PP_err_l2[1:p]))

  }

  p=6
  orth_b <- cbind(V[[1]],V[[2]],V[[3]],V[[4]],V[[5]],V[[6]])
  del <- (upper-lower)/T
  A <- del*Xkc%*%orth_b
  w <- ginv(t(A)%*%A)%*%t(A)%*%yc
  b_hat_pp[,p] <- w[1]*V[[1]]+w[2]*V[[2]]+w[3]*V[[3]]+w[4]*V[[4]]+w[5]*V[[5]]+w[6]*V[[6]]
  y_hat_pp[,p] <- w[1]*A[,1]+w[2]*A[,2]+w[3]*A[,3]+w[4]*A[,4]+w[5]*A[,5]+w[6]*A[,6]
  SSE_PP[p] <- sum((yc-y_hat_pp[,p])^2)/n
  AIC_PP[p] <- n*log(SSE_PP[p]) + 2*p
  BIC_PP[p] <- n*log(SSE_PP[p]) + p*log(n)
  y_pp_err_l2[p] <- sqrt(SSE_PP[p])

  if(!is.null(PPP_test)){
    yhat_test_PP[,p] <- y_bar + del*(X_k_test_c%*%b_hat_pp[,p])
    y_test_PP_err_l2[p] <- sqrt(sum((y_test-yhat_test_PP[,p])^2)/n)
  }

  if(p_pls == p){
    nbasis_PP <- order(BIC_PP[1:p])[1]
    return(list(b_hat = b_hat_pp[,1:p], y_hat = y_hat_pp[,1:p], BIC = BIC_PP[1:p], y_hat_test = yhat_test_PP[,1:p], basis = orth_b, rMSPE_test = y_test_PP_err_l2[1:p]))

  }


  p=7
  orth_b <- cbind(V[[1]],V[[2]],V[[3]],V[[4]],V[[5]],V[[6]],V[[7]])
  del <- (upper-lower)/T
  A <- del*Xkc%*%orth_b
  w <- ginv(t(A)%*%A)%*%t(A)%*%yc
  b_hat_pp[,p] <- w[1]*V[[1]]+w[2]*V[[2]]+w[3]*V[[3]]+w[4]*V[[4]]+w[5]*V[[5]]+w[6]*V[[6]]+w[7]*V[[7]]
  y_hat_pp[,p] <- w[1]*A[,1]+w[2]*A[,2]+w[3]*A[,3]+w[4]*A[,4]+w[5]*A[,5]+w[6]*A[,6]+w[7]*A[,7]
  SSE_PP[p] <- sum((yc-y_hat_pp[,p])^2)/n
  AIC_PP[p] <- n*log(SSE_PP[p]) + 2*p
  BIC_PP[p] <- n*log(SSE_PP[p]) + p*log(n)
  y_pp_err_l2[p] <- sqrt(SSE_PP[p])

  if(!is.null(PPP_test)){
    yhat_test_PP[,p] <- y_bar + del*(X_k_test_c%*%b_hat_pp[,p])
    y_test_PP_err_l2[p] <- sqrt(sum((y_test-yhat_test_PP[,p])^2)/n)
  }

  if(p_pls == p){
    nbasis_PP <- order(BIC_PP[1:p])[1]
    return(list(b_hat = b_hat_pp[,1:p], y_hat = y_hat_pp[,1:p], BIC = BIC_PP[1:p], y_hat_test = yhat_test_PP[,1:p], basis = orth_b, rMSPE_test = y_test_PP_err_l2[1:p]))

  }


  p=8
  orth_b <- cbind(V[[1]],V[[2]],V[[3]],V[[4]],V[[5]],V[[6]],V[[7]],V[[8]])
  del <- (upper-lower)/T
  A <- del*Xkc%*%orth_b
  w <- ginv(t(A)%*%A)%*%t(A)%*%yc
  b_hat_pp[,p] <- w[1]*V[[1]]+w[2]*V[[2]]+w[3]*V[[3]]+w[4]*V[[4]]+w[5]*V[[5]]+w[6]*V[[6]]+w[7]*V[[7]]+w[8]*V[[8]]
  y_hat_pp[,p] <- w[1]*A[,1]+w[2]*A[,2]+w[3]*A[,3]+w[4]*A[,4]+w[5]*A[,5]+w[6]*A[,6]+w[7]*A[,7]+w[8]*A[,8]
  SSE_PP[p] <- sum((yc-y_hat_pp[,p])^2)/n
  AIC_PP[p] <- n*log(SSE_PP[p]) + 2*p
  BIC_PP[p] <- n*log(SSE_PP[p]) + p*log(n)
  y_pp_err_l2[p] <- sqrt(SSE_PP[p])

  if(!is.null(PPP_test)){
    yhat_test_PP[,p] <- y_bar + del*(X_k_test_c%*%b_hat_pp[,p])
    y_test_PP_err_l2[p] <- sqrt(sum((y_test-yhat_test_PP[,p])^2)/n)
  }

  if(p_pls == p){
    nbasis_PP <- order(BIC_PP[1:p])[1]
    return(list(b_hat = b_hat_pp[,1:p], y_hat = y_hat_pp[,1:p], BIC = BIC_PP[1:p], y_hat_test = yhat_test_PP[,1:p], basis = orth_b, rMSPE_test = y_test_PP_err_l2[1:p]))

  }


  p=9
  orth_b <- cbind(V[[1]],V[[2]],V[[3]],V[[4]],V[[5]],V[[6]],V[[7]],V[[8]],V[[9]])
  del <- (upper-lower)/T
  A <- del*Xkc%*%orth_b
  w <- ginv(t(A)%*%A)%*%t(A)%*%yc
  b_hat_pp[,p] <- w[1]*V[[1]]+w[2]*V[[2]]+w[3]*V[[3]]+w[4]*V[[4]]+w[5]*V[[5]]+w[6]*V[[6]]+w[7]*V[[7]]+w[8]*V[[8]]+w[9]*V[[9]]
  y_hat_pp[,p] <- w[1]*A[,1]+w[2]*A[,2]+w[3]*A[,3]+w[4]*A[,4]+w[5]*A[,5]+w[6]*A[,6]+w[7]*A[,7]+w[8]*A[,8]+w[9]*A[,9]
  SSE_PP[p] <- sum((yc-y_hat_pp[,p])^2)/n
  AIC_PP[p] <- n*log(SSE_PP[p]) + 2*p
  BIC_PP[p] <- n*log(SSE_PP[p]) + p*log(n)
  y_pp_err_l2[p] <- sqrt(SSE_PP[p])

  if(!is.null(PPP_test)){
    yhat_test_PP[,p] <- y_bar + del*(X_k_test_c%*%b_hat_pp[,p])
    y_test_PP_err_l2[p] <- sqrt(sum((y_test-yhat_test_PP[,p])^2)/n)
  }

  if(p_pls == p){
    nbasis_PP <- order(BIC_PP[1:p])[1]
    return(list(b_hat = b_hat_pp[,1:p], y_hat = y_hat_pp[,1:p], BIC = BIC_PP[1:p], y_hat_test = yhat_test_PP[,1:p], basis = orth_b, rMSPE_test = y_test_PP_err_l2[1:p]))

  }


  p=10
  orth_b <- cbind(V[[1]],V[[2]],V[[3]],V[[4]],V[[5]],V[[6]],V[[7]],V[[8]],V[[9]],V[[10]])
  del <- (upper-lower)/T
  A <- del*Xkc%*%orth_b
  w <- ginv(t(A)%*%A)%*%t(A)%*%yc
  b_hat_pp[,p] <- w[1]*V[[1]]+w[2]*V[[2]]+w[3]*V[[3]]+w[4]*V[[4]]+w[5]*V[[5]]+w[6]*V[[6]]+w[7]*V[[7]]+w[8]*V[[8]]+w[9]*V[[9]]+w[10]*V[[10]]
  y_hat_pp[,p] <- w[1]*A[,1]+w[2]*A[,2]+w[3]*A[,3]+w[4]*A[,4]+w[5]*A[,5]+w[6]*A[,6]+w[7]*A[,7]+w[8]*A[,8]+w[9]*A[,9]+w[10]*A[,10]
  SSE_PP[p] <- sum((yc-y_hat_pp[,p])^2)/n
  AIC_PP[p] <- n*log(SSE_PP[p]) + 2*p
  BIC_PP[p] <- n*log(SSE_PP[p]) + p*log(n)
  y_pp_err_l2[p] <- sqrt(SSE_PP[p])

  if(!is.null(PPP_test)){
    yhat_test_PP[,p] <- y_bar + del*(X_k_test_c%*%b_hat_pp[,p])
    y_test_PP_err_l2[p] <- sqrt(sum((y_test-yhat_test_PP[,p])^2)/n)
  }


  nbasis_PP <- order(BIC_PP)[1]

  if(is.null(PPP_test)){
    yhat_test_PP <- NULL
    y_test_PP_err_l2 <- NULL
  }

  return(list(b_hat = b_hat_pp, y_hat = y_hat_pp, BIC = BIC_PP, y_hat_test = yhat_test_PP, basis = orth_b, rMSPE_test = y_test_PP_err_l2))
}




