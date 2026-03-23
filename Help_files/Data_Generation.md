# Data Generation
We suggest to call the following libraries before running the examples.

```r
library(spatstat)
library(fda)
library(nhppp)
library(mvtnorm)
library(Matrix)
library(pracma)
library(fdapace)
```

## Simulated Data `D_sim`

The R-code provided below generates the dataset `D_sim` in the R package. This corresponds to Case  of the simulation setting described in the manuscript. 


```r
set.seed(12345)

## Log intensity functions
n_basis <- 20
n_obs <- 100
n_test <- 100
bspl8 <- create.bspline.basis(rangeval = c(0,24), nbasis=n_basis)
eta <- 10
X_list_nobs <- lapply(1:(n_obs+n_test), function(i) fd(c(0, rnorm(1,12,4), rnorm(1,20,10), rnorm(1,20,10), rnorm(1,20,10),rnorm(1,20,10),rnorm(1,20,10),
                                                         rnorm(1,20,10), rnorm(1,20,10),rnorm(1,20,10), rnorm(1,20,10), rnorm(1,20,10), rnorm(1,20,10), rnorm(1,20,10), rnorm(1,20,10), rnorm(1,20,10), rnorm(1,20,10),rnorm(1,20,10), rnorm(1,20,10), rnorm(1,12,4))/eta + 2.8, bspl8))

X_list <- X_list_nobs[1:n_obs]
X_list_test <- X_list_nobs[(n_obs+1):(n_obs+n_test)]

## Coefficient function
b1 <-  fd(c(0, 1, 1, 1, 1,1,1,1,1, 1, -1, -1, -1,-1,-1,-1,-1,-1,-1,-1), bspl8)

## Response

### Training set
y <- rep(0,n_obs)
for(i in 1:n_obs){
  integrad <- function(t) predict.fd(b1, t)*predict.fd(X_list[[i]], t)
  y[i] <- integrate(integrad, lower = 0, upper = 24)$value + rnorm(1,0,1)
}

### Test set
y_test <- rep(0,n_test)
for(i in 1:n_test){
  integrad <- function(t) predict.fd(b1, t)*predict.fd(X_list_test[[i]], t)
  y_test[i] <- integrate(integrad, lower = 0, upper = 24)$value + rnorm(1,0,1)
}


n <- n_obs
T <- 100
lower <- 0
upper <- 24
del <- (upper-lower)/T
T_seq <- seq(0,24,length.out=T)
X <- matrix(unlist(lapply(X_list, function(x) predict.fd(x,T_seq))), nrow = n_obs, byrow = TRUE)
Xc <- t(t(X) - colMeans(X))
yc <- y - mean(y)
K_hat <- (1/n_obs)*t(Xc)%*%Xc

y_bar <- mean(y)

inprod_mat <- function(u1,u2,K_hat, lower, upper){
  del <- (upper-lower)/length(u1)
  value <- (del^2)*t(u1)%*%K_hat%*%u2
  return(value[1,1])
}

## Point process generation
PPP_obs_all <- lapply(X_list_nobs, function(x) draw_intensity(lambda = function(t) exp(predict.fd(x, t)), range_t = c(0,24), lambda_maj = c(intercept = 50000+max(exp(predict.fd(x, T_seq))))))
PPP_obs <- PPP_obs_all[1:n_obs]
PPP_test <- PPP_obs_all[(n_obs+1):(n_obs+n_test)]

## D_sim
b <- predict.fd(b1, T_seq)
X_list_all <- X_list_nobs
y_obs <- y

D_sim <- list("PPP_obs" = PPP_obs, "PPP_test" = PPP_test, "y_obs" =  y_obs, "y_test" = y_test, "b" = b, "X_list_all" =  X_list_all)

```
