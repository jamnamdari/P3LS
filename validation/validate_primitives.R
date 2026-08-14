## Validate that the new Rcpp primitives reproduce the original R
## implementations bit-for-bit (up to floating point tolerance).
## Run from the package root.

Rbin <- "C:/Program Files/R/R-4.6.0/bin/x64"
suppressMessages({
  library(spatstat)
  library(pracma)
  library(MASS)
})

## ---- 1. Load ORIGINAL implementations into their own environment ----
orig <- new.env()
sys.source("validation/P3LS_original.R", envir = orig)

## ---- 2. Load NEW (Rcpp-backed) package ----
suppressWarnings(suppressMessages(devtools::load_all(".", quiet = TRUE)))

set.seed(1)
cat("================ 1. edge_c ================\n")
cases <- list(
  c(s=0.5, h=1, lower=0, upper=24),
  c(s=0.1, h=1, lower=0, upper=24),   # near left boundary
  c(s=23.7, h=1, lower=0, upper=24),  # near right boundary
  c(s=12, h=2.5, lower=0, upper=24)
)
for (cc in cases) {
  o <- orig$edge_c(cc[["s"]], cc[["h"]], cc[["lower"]], cc[["upper"]])
  n <- edge_c_cpp(cc[["s"]], cc[["h"]], cc[["lower"]], cc[["upper"]])
  cat(sprintf("  s=%.2f h=%.2f -> orig=%.10f new=%.10f diff=%.3e\n", cc[["s"]], cc[["h"]], o, n, abs(o-n)))
}

cat("================ 2. lambda_fun (log intensity) ================\n")
load("data-raw/simulated_data.rda")
PP <- D_sim$PPP_obs[[1]]
T_seq <- seq(0, 24, length.out = 100)
h <- 1
lower <- 0; upper <- 24
orig_ll <- sapply(T_seq, function(s) log(orig$lambda_fun(PP, s, h)))
new_ll <- log_lambda_matrix_cpp(list(PP), T_seq, h, lower, upper)[1, ]
cat("  max abs diff:", max(abs(orig_ll - new_ll)), "\n")

cat("================ 3. inprod_mat ================\n")
Tn <- 40
u1 <- rnorm(Tn); u2 <- rnorm(Tn)
Kmat <- crossprod(matrix(rnorm(Tn*Tn), Tn, Tn))
o <- orig$inprod_mat(u1, u2, Kmat, 0, 24)
n <- inprod_mat_cpp(u1, u2, Kmat, 0, 24)
cat("  orig=", o, " new=", n, " diff=", abs(o-n), "\n")

cat("================ 4. Gram-Schmidt (GS) ================\n")
Ulist <- lapply(1:10, function(i) rnorm(Tn))
Vorig <- orig$GS(Ulist, Kmat, 0, 24)
Umat <- do.call(cbind, Ulist)
Vnew <- gram_schmidt_cpp(Umat, Kmat, 0, 24)
maxdiff <- max(sapply(1:10, function(j) max(abs(Vorig[[j]] - Vnew[, j]))))
cat("  max abs diff across 10 vectors:", maxdiff, "\n")
# orthonormality check (new)
G <- sapply(1:10, function(i) sapply(1:10, function(j) inprod_mat_cpp(Vnew[,i], Vnew[,j], Kmat, 0, 24)))
cat("  max |G - I| (orthonormality of new basis):", max(abs(G - diag(10))), "\n")

cat("================ 5. bin_counts ================\n")
nb <- 100
ends <- seq(0, 24, length.out = nb)
bins0 <- matrix(c(0, ends, ends, 0), nrow = nb+1, ncol = 2, byrow = FALSE)
bins <- bins0[-c(1, nb+1), ]
b_l <- bins[1,2]-bins[1,1]
PPP_sub <- D_sim$PPP_obs[1:10]
orig_counts <- t(sapply(PPP_sub, function(PP_ts) apply(bins, 1, function(t) sum(PP_ts>t[1] & PP_ts<=t[2]))))
new_counts <- bin_counts_matrix_cpp(PPP_sub, bins)
cat("  max abs diff:", max(abs(orig_counts - new_counts)), "\n")

cat("================ 6. Cov_estimator (fast path vs spatstat) ================\n")
## NOTE: sourcing the original file re-runs its top-level `library(nhppp)`
## call, and nhppp also exports a function called `ppp()` which then shadows
## spatstat::ppp() (a pre-existing conflict in the package's own library()
## calls, unrelated to this Rcpp change -- see chat notes). Detach it so the
## *original* fallback path can be exercised fairly for this comparison.
if ("package:nhppp" %in% search()) detach("package:nhppp", unload = TRUE, character.only = TRUE)

PPP_small <- D_sim$PPP_obs[1:15]
ngrid <- 40
bwd <- 1.5
t0 <- Sys.time()
Rorig <- orig$Cov_estimator(PPP_small, lbd=0, ubd=24, bwd=bwd, ngrid=ngrid, kern="epanechnikov")
t1 <- Sys.time()
Rnew <- Cov_estimator(PPP_small, lbd=0, ubd=24, bwd=bwd, ngrid=ngrid, kern="epanechnikov")
t2 <- Sys.time()
cat("  max abs diff R_X:", max(abs(Rorig - Rnew)), "\n")
cat("  orig time:", as.numeric(t1-t0, units="secs"), "s   new time:", as.numeric(t2-t1, units="secs"), "s\n")

cat("\nAll primitive validations complete.\n")
