\name{P3LS}
\alias{P3LS}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Point Process Partial Least Squares
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
Implements an extension of functional partial least square to point process data.
}
\usage{
P3LS(PPP_obs, PPP_test = NULL, y, y_test=NULL, h, p = 10, q, T, lbd, ubd, nb = 100)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{PPP_obs}{
%%     ~~Describe \code{x} here~~
A list containing vectors of tarining point process observed times.
}
  \item{PPP_test}{
%%     ~~Describe \code{x} here~~
A list containing vectors of testing point process observed times.
}
  \item{y}{
%%     ~~Describe \code{x} here~~
Training response values.
}
  \item{y_test}{
%%     ~~Describe \code{x} here~~
Testing response values.
}
  \item{h}{
%%     ~~Describe \code{x} here~~
Smoothing parameter in covariance function estimation.
}
  \item{p}{
%%     ~~Describe \code{x} here~~
Number of PLS basis used.
}
  \item{q}{
%%     ~~Describe \code{x} here~~
Number of basis functions used in log-intensity estimation.
}
  \item{T}{
%%     ~~Describe \code{x} here~~
Number of bins used for numerical integration.
}
  \item{lbd}{
%%     ~~Describe \code{x} here~~
Lower bound of the time interval on which the point process is observed.
}
  \item{ubd}{
%%     ~~Describe \code{x} here~~
Upper bound of the time interval on which the point process is observed.

}
  \item{nb}{
%%     ~~Describe \code{x} here~~
Number of bins in log-intensity estimation.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
\item{b_hat}{The object b_hat is a T by p matrix. Column k=1,...,p contains estimated coefficient function (using k PLS basis functions) evaluated on the grid of points of length T over [ldb,ubd].}
\item{y_hat}{The object y_hat is a T by p matrix. Column k=1,...,p contains predicted response values of the training data (using k PLS basis functions). }
\item{BIC}{The object BIC is a vector of length p. It`s k-th element, k=1,...,p, contains BIC critetian for the number of babsis functions (k) used in estimation of the coefficient function b.}
\item{y_hat_test}{{The object y_hat_test is a n_test by p matrix, where n_test is the number of testing subjects. Column k=1,...,p contains predicted response values of the testing data.}}
\item{basis}{The PLS basis functions. It computed the first p=10 basis functions evaluated on the grid of points of length T over [ldb,ubd]}
\item{rMSPE_test}{The object rMSPE_test is a vector of length p. It`s k-th element, k=1,...,p, contains the square root of the MSPE of the testing data set (using k PLS basis functions).}
}

\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
Jamshid Namdari
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
library(spatstat)
library(refund)
library(fda)
library(nhppp)
library(mvtnorm)
library(Matrix)
library(pracma)
library(fdapace)


#########################################
## P3LS
#########################################

Results_PP <- P3LS(D_sim$PPP_obs, D_sim$PPP_test, D_sim$y_obs, D_sim$y_test, h=2, p=10, q=15, T=100, lbd=0, ubd=24, nb = 100)

#########################################
## Plots
#########################################

## Plot of the estimated coefficient function
plot(seq(0,24,length.out=100), Results_PP$b_hat[,2], type='l', xlab = "t", ylab = "b(t)", main = "Coefficient function - Case 2")
points(seq(0,24,length.out=100), b, type="l", col="blue")
legend("topright", legend = c("estimated", "true"), lty=c(1,1), col = c("black", "blue"))


## Plot of Y_hat_Test vs y_test
plot(y_test, Results_PP$y_hat_test[,2], xlab="Response", ylab = "Estimated response", main = "Test data")

# rMSPE
Results_PP$rMSPE_test[2]

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
