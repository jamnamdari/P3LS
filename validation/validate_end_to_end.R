## End-to-end comparison: original vs Rcpp-optimized package, on D_sim.
suppressMessages({
  library(spatstat)
  library(pracma)
  library(MASS)
})

orig <- new.env()
sys.source("validation/P3LS_original.R", envir = orig)
if ("package:nhppp" %in% search()) detach("package:nhppp", unload = TRUE, character.only = TRUE)

suppressWarnings(suppressMessages(devtools::load_all(".", quiet = TRUE)))

load("data-raw/simulated_data.rda")

cmp <- function(a, b, label) {
  a <- as.numeric(unlist(a)); b <- as.numeric(unlist(b))
  d <- max(abs(a - b))
  cat(sprintf("  %-14s max abs diff = %.3e\n", label, d))
  invisible(d)
}

## ---------------- PLS_Kernel ----------------
cat("================ PLS_Kernel end-to-end (n=100, T=100) ================\n")
n_use <- 100
PPP_obs <- D_sim$PPP_obs[1:n_use]
PPP_test <- D_sim$PPP_test[1:n_use]
y <- D_sim$y_obs[1:n_use]
y_test <- D_sim$y_test[1:n_use]

t0 <- Sys.time()
r_orig <- orig$PLS_Kernel(PPP_obs, PPP_test, y, y_test, h=1, T=100, lbd=0, ubd=24)
t1 <- Sys.time()
r_new <- PLS_Kernel(PPP_obs, PPP_test, y, y_test, h=1, T=100, lbd=0, ubd=24)
t2 <- Sys.time()

cmp(r_orig$b_hat, r_new$b_hat, "b_hat")
cmp(r_orig$y_hat, r_new$y_hat, "y_hat")
cmp(r_orig$BIC, r_new$BIC, "BIC")
cat(sprintf("  orig time: %.2fs   new time: %.2fs   speedup: %.1fx\n",
            as.numeric(t1-t0,units="secs"), as.numeric(t2-t1,units="secs"),
            as.numeric(t1-t0,units="secs")/as.numeric(t2-t1,units="secs")))
## NOTE: y_hat_test / rMSPE_test are deliberately NOT compared to `orig`
## here. The original PLS_Kernel() used `t(Xkc_test) %*% b_hat_kr[,p]`
## (transposed) instead of `Xkc_test %*% b_hat_kr[,p]`, which only avoided
## erroring when n_test == T (as in this n=100,T=100 case) and even then
## computed the wrong quantity (summing over the wrong index). That has been
## fixed in the new code, so its y_hat_test/rMSPE_test legitimately differ
## from -- and are more correct than -- `orig`'s. See the n_test != T and
## before/after demonstrations below.
cat("  y_hat_test / rMSPE_test intentionally differ from `orig` -- see bugfix demo below\n")

cat("\n---- PLS_Kernel bugfix demo: y_hat_test with n_test != T ----\n")
n_small <- 20
r_bugfix <- PLS_Kernel(D_sim$PPP_obs[1:n_small], D_sim$PPP_test[1:n_small], D_sim$y_obs[1:n_small], D_sim$y_test[1:n_small], h=1, T=100, lbd=0, ubd=24)
r_orig_fail <- tryCatch({
  orig$PLS_Kernel(D_sim$PPP_obs[1:n_small], D_sim$PPP_test[1:n_small], D_sim$y_obs[1:n_small], D_sim$y_test[1:n_small], h=1, T=100, lbd=0, ubd=24)
  "no error (unexpected)"
}, error = function(e) paste("ERROR:", conditionMessage(e)))
cat("  orig$PLS_Kernel(n_test=20, T=100):", r_orig_fail, "\n")
cat("  new  PLS_Kernel(n_test=20, T=100): OK, dim(y_hat_test) =", paste(dim(r_bugfix$y_hat_test), collapse="x"), "\n")

## ---------------- P3LS (small subset: GLM loop untouched, kept slow) ----------------
## NOTE: P3LS()'s internal Cov_estimator() call hardcodes ngrid=100, and its
## GLM step does rbind(rep(1,nb), ef) then indexes X_k[,1:T] from the result
## -- both `nb` and `T` must equal that hardcoded 100 for the dimensions to
## conform (a pre-existing constraint of the original algorithm, not
## something introduced by this refactor -- see chat notes).
cat("\n================ P3LS end-to-end (n=30, T=100, nb=100, q=5, p=6) ================\n")
n_use2 <- 30
nb_use <- 100
PPP_obs2 <- D_sim$PPP_obs[1:n_use2]
PPP_test2 <- D_sim$PPP_test[1:n_use2]
y2 <- D_sim$y_obs[1:n_use2]
y_test2 <- D_sim$y_test[1:n_use2]

## devtools::load_all() above re-ran the package's own library(nhppp) call,
## which re-shadows spatstat::ppp() again (same pre-existing conflict noted
## earlier). Detach once more so the original reference path can run.
if ("package:nhppp" %in% search()) detach("package:nhppp", unload = TRUE, character.only = TRUE)

t0 <- Sys.time()
p_orig <- orig$P3LS(PPP_obs2, PPP_test2, y2, y_test2, h=1, p=6, q=5, T=100, lbd=0, ubd=24, nb=nb_use)
t1 <- Sys.time()
p_new <- P3LS(PPP_obs2, PPP_test2, y2, y_test2, h=1, p=6, q=5, T=100, lbd=0, ubd=24, nb=nb_use)
t2 <- Sys.time()

cmp(p_orig$b_hat, p_new$b_hat, "b_hat")
cmp(p_orig$y_hat, p_new$y_hat, "y_hat")
cmp(p_orig$BIC, p_new$BIC, "BIC")
cmp(p_orig$y_hat_test, p_new$y_hat_test, "y_hat_test")
cmp(p_orig$rMSPE_test, p_new$rMSPE_test, "rMSPE_test")
cat("  ncol(b_hat) orig/new:", ncol(p_orig$b_hat), "/", ncol(p_new$b_hat), "\n")
cat(sprintf("  orig time: %.2fs   new time: %.2fs   speedup: %.1fx\n",
            as.numeric(t1-t0,units="secs"), as.numeric(t2-t1,units="secs"),
            as.numeric(t1-t0,units="secs")/as.numeric(t2-t1,units="secs")))

## Check early-return quirk reproduced: p=1 and p=2 should both give 2 columns
cat("\n---- early-return quirk check (p_pls=1 vs p_pls=2) ----\n")
p1_new <- P3LS(PPP_obs2, PPP_test2, y2, y_test2, h=1, p=1, q=5, T=100, lbd=0, ubd=24, nb=nb_use)
p2_new <- P3LS(PPP_obs2, PPP_test2, y2, y_test2, h=1, p=2, q=5, T=100, lbd=0, ubd=24, nb=nb_use)
cat("  ncol(b_hat) for p_pls=1:", ncol(p1_new$b_hat), " (expect 2, matching original bug/behavior)\n")
cat("  ncol(b_hat) for p_pls=2:", ncol(p2_new$b_hat), " (expect 2)\n")

cat("\nAll end-to-end validations complete.\n")
