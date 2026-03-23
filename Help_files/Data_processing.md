# Data Processing

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

Please run the following code to get the objects `PPP_obs`, `PPP_test`, `y`, and `y_test`. 


```r
n_obs0 <- length(D)
n_obs <- 100
n_test <- n_obs0 - n_obs


count2pp_temp <- function(D1){
  PP1 <- c()
  ub <- 0
  for(k in 1:nrow(D1)){
    ub <- ifelse(k == nrow(D1),2700,D1[k+1,2][[1]])
    PP1 <- c(PP1, runif(D1[k,3][[1]],D1[k,2][[1]],ub)/60)
  }
  return(PP1)
}

PPP_obs_all <- lapply(D[1:n_obs0], count2pp_temp)
PPP_obs <- PPP_obs_all[1:n_obs]
PPP_test <- PPP_obs_all[(n_obs+1):n_obs0]

y_all <- sapply(D, function(x) x$AT[1])
y <- y_all[1:n_obs]
y_test <- y_all[(n_obs+1):n_obs0]

```
