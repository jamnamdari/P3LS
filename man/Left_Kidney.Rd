\name{Left_Kidney}
\alias{Left_Kidney}
\docType{data}
\title{
%%   ~~ data name/kind ... ~~
Left Kidney Renogram Data
}
\description{
%%  ~~ A concise (1-5 lines) description of the dataset. ~~
This dataset contains the Renogram data and an experts ratings of obstruction status of the left kidney of patients who were referred to the clinic suspected with kidney obstruction.
}
\usage{data("Left_Kidney")}
\format{
  A list of data frames with 127 observations on the following 4 variables.
  \describe{
    \item{\code{Patient.Code.Name}}{a character vector.}
    \item{\code{Time_Interval_Stamp}}{is an integer vector contains end point of imaging time intevals.}
    \item{\code{renogram_value}}{is an integer vector containing photone counts.}
    \item{\code{AT}}{is a numeric vector containing an expert`s ratings ranges from -1.0 to 1.0 with higher scores indicating higher likelihood of obstruction.}

  }
}
\details{
%%  ~~ If necessary, more details than the __description__ above ~~
}
\source{
%%  ~~ reference to a publication or URL from which the data were obtained ~~
}
\references{
%%  ~~ possibly secondary sources and usages ~~
}
\examples{
library(mvtnorm)
library(Matrix)
library(fda)
library(nhppp)
library(spatstat)
library(refund)
library(dplyr)
library(pracma)
library(fdapace)

## Loading the left kidney data.
data("Left_Kidney")

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

####################################################
## P3LS
####################################################
Results_Left <- P3LS:::P3LS(PPP_obs, PPP_test, y, y_test, h=2, p=10, q=10, T=100, lbd=0, ubd=45, nb = 100)

####################################################
## Plots
####################################################

# Plotof the estimated coefficient function
plot(Results_Left$b_hat[,2], lty=1, main = "Coefficient Function - Left Kidney", type="l", ylab="", xlab = "time index", ylim=c(-.02,.025))

# Plot of the first PLS basis function
plot(Results_Left$basis[,1], xlab=", ylab=", main="Frist PLS basis function - Left Kidnay")

# Plot of the second PLS basis function
plot(Results_Left$basis[,2], xlab=", ylab=", main="Second PLS basis function - Left Kidnay")

# Plot of the square root of MSPE for p=1,...,10
plot(1:10, Results_Left$rMSPE_test, xlab="Number of PLS basis functions - Left Kidney", ylab="rMSPE", main="Square root of MSPE")

}
\keyword{datasets}
