## Generates all data (timing, correctness diffs, applied-validity results)
## for the comparison report. Writes validation/report_data.json.
suppressMessages({
  library(spatstat); library(pracma); library(MASS); library(jsonlite)
})

orig <- new.env()
sys.source("validation/P3LS_original.R", envir = orig)
if ("package:nhppp" %in% search()) detach("package:nhppp", unload = TRUE, character.only = TRUE)

suppressWarnings(suppressMessages(devtools::load_all(".", quiet = TRUE)))
if ("package:nhppp" %in% search()) detach("package:nhppp", unload = TRUE, character.only = TRUE)

load("data-raw/simulated_data.rda")

report <- list()
timeit <- function(expr) {
  t0 <- Sys.time()
  val <- force(expr)
  list(value = val, secs = as.numeric(Sys.time() - t0, units = "secs"))
}

## ================= Section A: primitive-level correctness =================
cat("Section A: primitives\n")
A <- list()

cases <- list(c(0.5,1,0,24), c(0.1,1,0,24), c(23.7,1,0,24), c(12,2.5,0,24))
A$edge_c <- lapply(cases, function(cc) {
  o <- orig$edge_c(cc[1], cc[2], cc[3], cc[4])
  n <- edge_c_cpp(cc[1], cc[2], cc[3], cc[4])
  list(s=cc[1], h=cc[2], orig=o, new=n, diff=abs(o-n))
})

PP <- D_sim$PPP_obs[[1]]
T_seq100 <- seq(0, 24, length.out = 100)
lower <- 0; upper <- 24  # orig$lambda_fun()'s internal edge_c() call relies on these as free variables
orig_ll <- sapply(T_seq100, function(s) log(orig$lambda_fun(PP, s, 1)))
new_ll <- log_lambda_matrix_cpp(list(PP), T_seq100, 1, 0, 24)[1, ]
A$lambda_fun_maxdiff <- max(abs(orig_ll - new_ll))

set.seed(1)
Tn <- 40
u1 <- rnorm(Tn); u2 <- rnorm(Tn)
Kmat <- crossprod(matrix(rnorm(Tn*Tn), Tn, Tn))
A$inprod_diff <- abs(orig$inprod_mat(u1,u2,Kmat,0,24) - inprod_mat_cpp(u1,u2,Kmat,0,24))

Ulist <- lapply(1:10, function(i) rnorm(Tn))
Vorig <- orig$GS(Ulist, Kmat, 0, 24)
Vnew <- gram_schmidt_cpp(do.call(cbind, Ulist), Kmat, 0, 24)
A$gs_maxdiff <- max(sapply(1:10, function(j) max(abs(Vorig[[j]] - Vnew[, j]))))

nb <- 100
ends <- seq(0, 24, length.out = nb)
bins0 <- matrix(c(0, ends, ends, 0), nrow = nb+1, ncol = 2, byrow = FALSE)
bins <- bins0[-c(1, nb+1), ]
PPP_sub <- D_sim$PPP_obs[1:10]
orig_counts <- t(sapply(PPP_sub, function(PP_ts) apply(bins, 1, function(t) sum(PP_ts>t[1] & PP_ts<=t[2]))))
new_counts <- bin_counts_matrix_cpp(PPP_sub, bins)
A$bin_counts_maxdiff <- max(abs(orig_counts - new_counts))

cov_orig <- timeit(orig$Cov_estimator(D_sim$PPP_obs[1:15], lbd=0, ubd=24, bwd=1.5, ngrid=40, kern="epanechnikov"))
cov_new  <- timeit(Cov_estimator(D_sim$PPP_obs[1:15], lbd=0, ubd=24, bwd=1.5, ngrid=40, kern="epanechnikov"))
A$cov_estimator_maxdiff <- max(abs(cov_orig$value - cov_new$value))
A$cov_estimator_time <- list(orig = cov_orig$secs, new = cov_new$secs)

report$primitives <- A

## ================= Section B: PLS_Kernel timing scaling =================
## `orig`'s y_hat_kr is allocated matrix(0, nrow=T, ncol=dd) instead of
## nrow=n_obs (the third bug found while building this report -- see chat
## notes): it errors outright whenever n_obs does not evenly divide T, and
## SILENTLY RECYCLES (i.e. produces wrong, uncaught results) whenever it
## does, e.g. n_obs=20 with T=100. It is only actually valid when
## n_obs == T. So the timing SCALING curve below uses `new` only (the one
## implementation that is correct at every n); the orig-vs-new comparison is
## reported separately at the one n where orig is valid (n=100), plus an
## explicit demonstration of both failure modes (hard error at n=40, silent
## corruption at n=20).
cat("Section B: PLS_Kernel scaling\n")
ns <- c(20, 40, 60, 80, 100)
B <- list(new_scaling = list())
for (n_use in ns) {
  PPo <- D_sim$PPP_obs[1:n_use]
  yo <- D_sim$y_obs[1:n_use]
  rn <- timeit(PLS_Kernel(PPo, NULL, yo, NULL, h=1, T=100, lbd=0, ubd=24))
  B$new_scaling[[as.character(n_use)]] <- list(n = n_use, new_secs = rn$secs)
  cat(sprintf("  n=%d  new=%.4fs\n", n_use, rn$secs))
}

## Valid comparison point: n_obs == T == 100 (the only config where orig's
## y_hat_kr allocation isn't broken).
PPo <- D_sim$PPP_obs; yo <- D_sim$y_obs
ro100 <- timeit(orig$PLS_Kernel(PPo, NULL, yo, NULL, h=1, T=100, lbd=0, ubd=24))
rn100 <- timeit(PLS_Kernel(PPo, NULL, yo, NULL, h=1, T=100, lbd=0, ubd=24))
B$valid_comparison_n100 <- list(
  orig_secs = ro100$secs, new_secs = rn100$secs,
  b_hat_maxdiff = max(abs(ro100$value$b_hat - rn100$value$b_hat)),
  y_hat_maxdiff = max(abs(ro100$value$y_hat - rn100$value$y_hat)),
  BIC_maxdiff = max(abs(ro100$value$BIC - rn100$value$BIC))
)
cat(sprintf("  [valid n=100] orig=%.3fs new=%.4fs speedup=%.1fx\n", ro100$secs, rn100$secs, ro100$secs/rn100$secs))

## Failure-mode demonstrations
PPo20 <- D_sim$PPP_obs[1:20]; yo20 <- D_sim$y_obs[1:20]
ro20 <- orig$PLS_Kernel(PPo20, NULL, yo20, NULL, h=1, T=100, lbd=0, ubd=24)
rn20 <- PLS_Kernel(PPo20, NULL, yo20, NULL, h=1, T=100, lbd=0, ubd=24)
B$silent_corruption_n20 <- list(
  n = 20, orig_ran_without_error = TRUE,
  orig_BIC = as.numeric(ro20$BIC), new_BIC = as.numeric(rn20$BIC),
  note = "orig does not error (100 %% 20 == 0, so R recycles y_hat_kr's 20-element column to fill 100 slots) but SSE/AIC/BIC are computed over the recycled vector against a length-20 yc, silently producing incoherent numbers."
)
n40_err <- tryCatch({
  orig$PLS_Kernel(D_sim$PPP_obs[1:40], NULL, D_sim$y_obs[1:40], NULL, h=1, T=100, lbd=0, ubd=24)
  NA
}, error = function(e) conditionMessage(e))
B$hard_error_n40 <- list(n = 40, orig_error = n40_err)
cat("  [n=20] orig BIC (silently corrupted):", paste(round(ro20$BIC,2), collapse=", "), "\n")
cat("  [n=20] new  BIC (correct):           ", paste(round(rn20$BIC,2), collapse=", "), "\n")
cat("  [n=40] orig error:", n40_err, "\n")

report$pls_kernel_scaling <- B

## ================= Section C: Cov_estimator timing scaling =================
cat("Section C: Cov_estimator scaling\n")
Cres <- list(by_ngrid = list(), by_n = list())
ngrids <- c(20, 40, 60, 80, 100)
for (ng in ngrids) {
  ro <- timeit(orig$Cov_estimator(D_sim$PPP_obs[1:15], lbd=0, ubd=24, bwd=1.5, ngrid=ng, kern="epanechnikov"))
  rn <- timeit(Cov_estimator(D_sim$PPP_obs[1:15], lbd=0, ubd=24, bwd=1.5, ngrid=ng, kern="epanechnikov"))
  Cres$by_ngrid[[as.character(ng)]] <- list(ngrid = ng, orig_secs = ro$secs, new_secs = rn$secs,
                                             maxdiff = max(abs(ro$value - rn$value)))
  cat(sprintf("  ngrid=%d  orig=%.3fs new=%.3fs\n", ng, ro$secs, rn$secs))
}
nsubj <- c(5, 10, 15, 20, 25)
for (ns2 in nsubj) {
  ro <- timeit(orig$Cov_estimator(D_sim$PPP_obs[1:ns2], lbd=0, ubd=24, bwd=1.5, ngrid=50, kern="epanechnikov"))
  rn <- timeit(Cov_estimator(D_sim$PPP_obs[1:ns2], lbd=0, ubd=24, bwd=1.5, ngrid=50, kern="epanechnikov"))
  Cres$by_n[[as.character(ns2)]] <- list(n = ns2, orig_secs = ro$secs, new_secs = rn$secs,
                                          maxdiff = max(abs(ro$value - rn$value)))
  cat(sprintf("  n=%d  orig=%.3fs new=%.3fs\n", ns2, ro$secs, rn$secs))
}
report$cov_estimator_scaling <- Cres

## ================= Section D: P3LS timing scaling =================
cat("Section D: P3LS scaling\n")
D <- list()
ns_p3ls <- c(10, 20, 30, 40, 50)
for (n_use in ns_p3ls) {
  PPo <- D_sim$PPP_obs[1:n_use]; PPt <- D_sim$PPP_test[1:n_use]
  yo <- D_sim$y_obs[1:n_use]; yt <- D_sim$y_test[1:n_use]
  ro <- timeit(suppressWarnings(orig$P3LS(PPo, PPt, yo, yt, h=1, p=6, q=5, T=100, lbd=0, ubd=24, nb=100)))
  rn <- timeit(suppressWarnings(P3LS(PPo, PPt, yo, yt, h=1, p=6, q=5, T=100, lbd=0, ubd=24, nb=100)))
  D[[as.character(n_use)]] <- list(
    n = n_use, orig_secs = ro$secs, new_secs = rn$secs,
    b_hat_maxdiff = max(abs(ro$value$b_hat - rn$value$b_hat)),
    rMSPE_maxdiff = max(abs(ro$value$rMSPE_test - rn$value$rMSPE_test))
  )
  cat(sprintf("  n=%d  orig=%.3fs new=%.3fs speedup=%.1fx\n", n_use, ro$secs, rn$secs, ro$secs/rn$secs))
}
report$p3ls_scaling <- D

## ================= Section E: applied validity (full n=100) =================
cat("Section E: applied validity\n")
r_full <- PLS_Kernel(D_sim$PPP_obs, D_sim$PPP_test, D_sim$y_obs, D_sim$y_test, h=1, T=100, lbd=0, ubd=24)
best_p <- order(r_full$BIC)[1]
E <- list(
  T_seq = T_seq100,
  b_true = as.numeric(D_sim$b),
  b_hat_best = as.numeric(r_full$b_hat[, best_p]),
  best_p = best_p,
  BIC = as.numeric(r_full$BIC),
  y_test = as.numeric(D_sim$y_test),
  y_hat_test_best = as.numeric(r_full$y_hat_test[, best_p]),
  rMSPE_test_new = as.numeric(r_full$rMSPE_test)
)

## "before" (buggy) rMSPE-by-p curve for the same n=100,T=100 case (this is
## the one configuration where orig's PLS_Kernel doesn't error, since
## n_test==T==100 there, but it still uses the wrong formula)
r_full_orig <- orig$PLS_Kernel(D_sim$PPP_obs, D_sim$PPP_test, D_sim$y_obs, D_sim$y_test, h=1, T=100, lbd=0, ubd=24)
E$rMSPE_test_orig_buggy <- as.numeric(r_full_orig$rMSPE_test)

## demonstrate the crash on non-square n_test != T, which orig cannot do at all
demo_ok <- tryCatch({
  orig$PLS_Kernel(D_sim$PPP_obs[1:20], D_sim$PPP_test[1:20], D_sim$y_obs[1:20], D_sim$y_test[1:20], h=1, T=100, lbd=0, ubd=24)
  TRUE
}, error = function(e) FALSE)
E$orig_fails_when_ntest_neq_T <- !demo_ok
new_demo <- PLS_Kernel(D_sim$PPP_obs[1:20], D_sim$PPP_test[1:20], D_sim$y_obs[1:20], D_sim$y_test[1:20], h=1, T=100, lbd=0, ubd=24)
E$new_works_when_ntest_neq_T <- all(dim(new_demo$y_hat_test) == c(20, 10))

report$applied_validity <- E

## ================= Section F: nhppp/spatstat shadowing demo =================
cat("Section F: shadowing bug demo\n")
Fsec <- list()
Fsec$nhppp_attached <- "package:nhppp" %in% search()
r_gauss <- tryCatch({
  Cov_estimator(D_sim$PPP_obs[1:8], lbd=0, ubd=24, bwd=1.5, ngrid=20, kern="gaussian")
  TRUE
}, error = function(e) FALSE)
Fsec$new_fallback_kernel_ok_under_shadowing <- r_gauss
report$shadowing_demo <- Fsec

writeLines(toJSON(report, auto_unbox = TRUE, digits = 10, na = "null"), "validation/report_data.json")
cat("\nWrote validation/report_data.json\n")
