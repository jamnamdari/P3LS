library(RSpectra)
library(mvtnorm)
library(Matrix)
library(lattice)
library(fda)
library(nhppp)
library(pracma)
library(ggplot2)
#library(mnormt)

#library(R.matlab)
library(spatstat)
#library(compiler)
#library(parallel)
#library(foreach)
#library(doParallel)
#library(doSNOW)
#library(refund)
library(fdapace)

##################################################################################
##
## NOTE ON THIS FILE (Rcpp optimization pass)
## ------------------------------------------
## The numerically heavy inner loops of this package have been moved to C++
## (see src/p3ls.cpp / R/RcppExports.R):
##   - inprod_mat() / GS()      -> inprod_mat_cpp() / gram_schmidt_cpp()
##   - lambda_fun()/Ker()/edge_c() grid evaluation in PLS_Kernel()
##                               -> log_lambda_matrix_cpp()
##   - bin counting in P3LS()   -> bin_counts_matrix_cpp()
##   - the spatstat::crosspairs-based inner loop of Cov_estimator()
##                               -> cov_estimator_core_cpp() (kernel="epanechnikov" only;
##                                  any other kernel keeps using the original
##                                  spatstat-based implementation as a fallback)
## The R-level function signatures and return values are unchanged; the
## copy-pasted "p = 1, 2, 3, ..." blocks in PLS_Kernel() and P3LS() have been
## collapsed into ordinary for() loops (generalizing GS() to an arbitrary
## number of vectors, rather than hardcoding 10), but the control flow
## (including the p_pls<=2 / p_pls==p early-return checks in P3LS(), and the
## inv() vs ginv() threshold at p<=4 vs p>=5) is preserved exactly.
##
##################################################################################

inprod_mat <- function(u1,u2,K_hat, lower, upper){
  inprod_mat_cpp(as.numeric(u1), as.numeric(u2), as.matrix(K_hat), lower, upper)
}



GS <- function(u_list,K_hat, lower, upper){
  Umat <- do.call(cbind, u_list)
  Vmat <- gram_schmidt_cpp(Umat, as.matrix(K_hat), lower, upper)
  lapply(seq_len(ncol(Vmat)), function(j) Vmat[, j])
}


### **Covariance estimation by point process method**

Cov_estimator <- function(PROCESS, lbd, ubd, bwd, ngrid, kern="epanechnikov"){

  lower <- lbd
  upper <- ubd
  n <- length(PROCESS)
  grids <- seq(lbd,ubd,l=ngrid)

  if (kern == "epanechnikov") {
    ## Fast path: closed-form Epanechnikov kernel, no spatstat dependency.
    ## Edge correction reproduces spatstat::pkernel(..., kernel="epanechnikov")
    ## exactly, including the original code's asymmetric use of `sd` between
    ## the lower- and upper-boundary terms (sd = 1 on the lower term, sd =
    ## sqrt(1/5) on the upper term) -- preserved for numerical equivalence.
    pkernel_epan <- function(x, sd){
      a <- sd*sqrt(5)
      y <- x/a
      p <- ifelse(y < -1, 0, ifelse(y > 1, 1, (2+3*y-y^3)/4))
      p
    }
    edge <- pkernel_epan((grids-lbd)/bwd, sd=1) - pkernel_epan((grids-ubd)/bwd, sd=sqrt(1/5))
    edge <- outer(edge,edge,FUN="*")

    core <- cov_estimator_core_cpp(PROCESS, grids, bwd)
    A_sum <- core$A_sum
    Gsum <- core$Gsum
    tmp1_sum <- as.numeric(core$tmp1_sum)

    A2 <- A_sum + Gsum
    C_pooled <- outer(tmp1_sum, tmp1_sum)

    A <- A_sum/(n*edge)
    C <- (C_pooled - A2)/((n*(n-1))*edge)

    R_X <- log(A/C)
    return(R_X)
  }

  ## ---- Fallback: original spatstat-based implementation for any other kernel ----
  ## NOTE: ppp()/owin()/crosspairs()/dkernel()/pkernel() are called via
  ## explicit spatstat.geom::/spatstat.univar:: namespace qualification
  ## rather than as bare names. nhppp (loaded via library(nhppp) at the top
  ## of this file) also exports a function called ppp() with a different
  ## signature; depending on package attach order it can shadow
  ## spatstat.geom::ppp() on the search path and make this fallback path
  ## fail with "unused arguments". Qualifying the calls removes the
  ## ambiguity regardless of load order.
  edge <- spatstat.univar::pkernel((grids-lbd)/bwd,kernel = kern)-spatstat.univar::pkernel((grids-ubd)/bwd,kernel = kern,sd=sqrt(1/5))
  edge <- outer(edge,edge,FUN="*")

  Kh <- function(t){
    spatstat.univar::dkernel(t/bwd,kernel = kern,sd=sqrt(1/5))/bwd
  }
  A2 <- A <- B <- C <- D <- matrix(0,ngrid,ngrid)

  m <- 1
  A2.out <- list()
  A2.out.y <- vector(mode = "list", length = m)
  for(i in 1:n){
    tmp.A2 <- matrix(0,ngrid,ngrid)
    for(j in 1:m){
      process <- PROCESS[[i]]
      if(!is.null(process)){
        cp <- spatstat.geom::crosspairs(spatstat.geom::ppp(x=process,y=process*0, window = spatstat.geom::owin(c(lbd,ubd),c(lbd,ubd))),spatstat.geom::ppp(x=grids,y=grids*0, window = spatstat.geom::owin(c(lbd,ubd),c(lbd,ubd))),rmax = bwd,what = "ijd")
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
    cp <- spatstat.geom::crosspairs(spatstat.geom::ppp(x=process,y=process*0, window = spatstat.geom::owin(c(lbd,ubd),c(lbd,ubd))),spatstat.geom::ppp(x=grids,y=grids*0, window = spatstat.geom::owin(c(lbd,ubd),c(lbd,ubd))),rmax = bwd,what = "ijd")
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
  edge_c_cpp(s, h, lower, upper)
}

#' @export
lambda_fun <- function(PP, s,h){
  sum(Ker(s,PP,h))/edge_c(s,h, lower, upper)
}


#' @export
PLS_Kernel <- function(PPP_obs, PPP_test = NULL, y, y_test=NULL, h, T, lbd, ubd){
  lower <- lbd
  upper <- ubd

  n_obs <- length(PPP_obs)
  if(!is.null(PPP_test)){
    n_test <- length(PPP_test)
  }
  n <- n_obs
  dd <- 10
  b_hat_kr  <- matrix(0,nrow = T, ncol=dd)
  ## y_hat_kr holds one predicted value per TRAINING SUBJECT (n_obs rows),
  ## the same T/n mixup as yhat_test_K below (and just as silent: for
  ## n_obs values that happen to divide T exactly, R recycles instead of
  ## erroring, corrupting SSE_K/AIC_K/BIC_K without a warning). Matches the
  ## correct, analogous `y_hat_pp <- matrix(0,nrow=n_obs,...)` in P3LS().
  y_hat_kr  <- matrix(0,nrow = n_obs, ncol=dd)
  y_kr_err_l2  <- rep(0,dd)
  ## yhat_test_K holds one predicted value per TEST SUBJECT (n_test rows),
  ## not per grid point -- the original code allocated this as `nrow = T`,
  ## which only happened to work when n_test == T by coincidence (as in the
  ## packaged example, where both are 100) and threw "non-conformable
  ## arguments" otherwise. When PPP_test is NULL, n_test is never assigned,
  ## so fall back to T as an (unused, discarded below) placeholder.
  yhat_test_K  <- matrix(0,nrow = if(is.null(PPP_test)) T else n_test, ncol=dd)
  y_test_K_err_l2  <- rep(0,dd)

  yc <- y - mean(y)
  y_bar <- mean(y)


  del <- (upper-lower)/T
  T_seq <- seq(lbd,ubd,length.out=T)

  X_k <- log_lambda_matrix_cpp(PPP_obs, T_seq, h, lower, upper)
  Xkc <- t(t(X_k) - colMeans(X_k))

  Kk_hat <- (1/n_obs)*t(Xkc)%*%Xkc
  K_b <- (1/n_obs)*t(Xkc)%*%yc

  ## Build K_b, K_hat %*% K_b, K_hat^2 %*% K_b, ..., K_hat^9 %*% K_b
  Kpow_b <- vector("list", dd)
  Kpow_b[[1]] <- K_b
  for(k in 2:dd){
    Kpow_b[[k]] <- (del)*Kk_hat%*%Kpow_b[[k-1]]
  }

  U <- Kpow_b
  V <- GS(U,Kk_hat, lower, upper)
  V_mat <- do.call(cbind, V)

  AIC_K <- rep(0,dd)
  BIC_K <- rep(0,dd)
  SSE_K <- rep(0,dd)

  if(!is.null(PPP_test)){
    X_k_test <- log_lambda_matrix_cpp(PPP_test, T_seq, h, lower, upper)
    Xkc_test <- t(t(X_k_test) - colMeans(X_k_test))
  }

  del <- (upper-lower)/T
  orth_b <- V_mat
  for(p in 1:dd){
    orth_b <- V_mat[,1:p,drop=FALSE]
    A <- del*Xkc%*%orth_b
    if(p <= 4){
      w <- inv(t(A)%*%A)%*%t(A)%*%yc
    } else {
      w <- ginv(t(A)%*%A)%*%t(A)%*%yc
    }
    b_hat_kr[,p] <- orth_b%*%w
    y_hat_kr[,p] <- A%*%w
    SSE_K[p] <- sum((yc-y_hat_kr[,p])^2)/n
    AIC_K[p] <- n*log(SSE_K[p]) + 2*p
    BIC_K[p] <- n*log(SSE_K[p]) + p*log(n)
    y_kr_err_l2[p] <- sqrt(SSE_K[p])

    if(!is.null(PPP_test)){
      ## NOT transposed: Xkc_test is n_test x T and b_hat_kr[,p] is a
      ## length-T coefficient function, so Xkc_test %*% b_hat_kr[,p] gives
      ## one prediction per test subject -- matching the analogous (correct)
      ## line in P3LS(). The original `t(Xkc_test) %*% b_hat_kr[,p]` summed
      ## over the wrong index (treating b_hat_kr's grid-point index as a
      ## subject index) and only avoided erroring outright when T == n_test
      ## happened to hold.
      yhat_test_K[,p] <- y_bar + del*(Xkc_test%*%b_hat_kr[,p])
      y_test_K_err_l2[p] <- sqrt(sum((y_test-yhat_test_K[,p])^2)/n)
    }
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
  library(pracma)
  library(MASS)
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

  ## Number of data points (events) in each bin, for every subject at once
  bin_counts_obs <- bin_counts_matrix_cpp(PPP_obs, bins)/b_l

  X_k <- matrix(0,nrow = n_obs, ncol = T)
  for(ell in 1:n_obs){
    Dm <- data.frame(counts = bin_counts_obs[ell,], phi = PHIsm)
    Pm_model <- glm(counts~.+1  , family = poisson(link = "log"), data = Dm)
    X_k[ell,] <- colSums(diag(Pm_model$coefficients)%*%rbind(rep(1,nb),ef))
  }
  Xkc <- t(t(X_k) - colMeans(X_k))

  # test set
  if(!is.null(PPP_test)){
    bin_counts_test <- bin_counts_matrix_cpp(PPP_test, bins)/b_l
    X_k_test <- matrix(0,nrow = n_test, ncol = T)
    for(ell in 1:n_test){
      Dm <- data.frame(counts = bin_counts_test[ell,], phi = PHIsm)
      Pm_model <- glm(counts~.+1  , family = poisson(link = "log"), data = Dm)
      X_k_test[ell,] <- colSums(diag(Pm_model$coefficients)%*%rbind(rep(1,nb),ef))
    }
    X_k_test_c <- t(t(X_k_test) - colMeans(X_k_test))
  }

  K_hat_e <- K_hat_pp1_e
  K_b <- (1/n_obs)*t(Xkc)%*%yc

  Kpow_b <- vector("list", dd)
  Kpow_b[[1]] <- K_b
  for(k in 2:dd){
    Kpow_b[[k]] <- (del)*K_X%*%Kpow_b[[k-1]]
  }

  U <- Kpow_b
  V <- GS(U,K_X, lower, upper)
  V_mat <- do.call(cbind, V)

  del <- (upper-lower)/T
  orth_b <- V_mat
  for(p in 1:dd){
    orth_b <- V_mat[,1:p,drop=FALSE]
    A <- del*Xkc%*%orth_b
    if(p <= 4){
      w <- inv(t(A)%*%A)%*%t(A)%*%yc
    } else {
      w <- ginv(t(A)%*%A)%*%t(A)%*%yc
    }
    b_hat_pp[,p] <- orth_b%*%w
    y_hat_pp[,p] <- A%*%w
    SSE_PP[p] <- sum((yc-y_hat_pp[,p])^2)/n
    AIC_PP[p] <- n*log(SSE_PP[p]) + 2*p
    BIC_PP[p] <- n*log(SSE_PP[p]) + p*log(n)
    y_pp_err_l2[p] <- sqrt(SSE_PP[p])

    if(!is.null(PPP_test)){
      yhat_test_PP[,p] <- y_bar + del*(X_k_test_c%*%b_hat_pp[,p])
      y_test_PP_err_l2[p] <- sqrt(sum((y_test-yhat_test_PP[,p])^2)/n)
    }

    ## Early-return behaviour reproduced exactly from the original
    ## implementation: p_pls %in% c(1,2) both stop right after the p=2
    ## block (so both return 2 columns); p_pls %in% 3:9 stop right after
    ## their own block; p_pls >= 10 falls through to the unconditional
    ## return after the full loop.
    if(p == 2 && p_pls <= 2){
      nbasis_PP <- order(BIC_PP[1:p])[1]
      return(list(b_hat = b_hat_pp[,1:p], y_hat = y_hat_pp[,1:p], BIC = BIC_PP[1:p], y_hat_test = yhat_test_PP[,1:p], basis = orth_b, rMSPE_test = y_test_PP_err_l2[1:p]))
    }
    if(p >= 3 && p <= 9 && p_pls == p){
      nbasis_PP <- order(BIC_PP[1:p])[1]
      return(list(b_hat = b_hat_pp[,1:p], y_hat = y_hat_pp[,1:p], BIC = BIC_PP[1:p], y_hat_test = yhat_test_PP[,1:p], basis = orth_b, rMSPE_test = y_test_PP_err_l2[1:p]))
    }
  }

  nbasis_PP <- order(BIC_PP)[1]

  if(is.null(PPP_test)){
    yhat_test_PP <- NULL
    y_test_PP_err_l2 <- NULL
  }

  return(list(b_hat = b_hat_pp, y_hat = y_hat_pp, BIC = BIC_PP, y_hat_test = yhat_test_PP, basis = orth_b, rMSPE_test = y_test_PP_err_l2))
}
