// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <algorithm>
#include <vector>

using namespace Rcpp;

// ============================================================================
// Epanechnikov kernel helpers (closed form)
//
// The R-level Ker()/edge_c() pair (used by lambda_fun() and PLS_Kernel()) is
// a plain Epanechnikov kernel with bandwidth h:
//     Ker(s,u,h) = (1/h)*(3/4)*(1-((s-u)/h)^2)   for |s-u| < h,  else 0
// edge_c(s,h,lower,upper) = integral of Ker(s, . ,h) over [lower,upper].
// That integral has a closed form (antiderivative of the Epanechnikov
// density), so it does not need numerical quadrature (integrate()).
// ============================================================================

// Antiderivative of the standard Epanechnikov density (3/4)(1-v^2) on
// |v|<1, clipped outside: F(v) in [-1/2, 1/2].
static inline double epan_F(double v) {
  if (v <= -1.0) return -0.5;
  if (v >= 1.0) return 0.5;
  return (3.0 * v - v * v * v) / 4.0;
}

// [[Rcpp::export]]
double edge_c_cpp(double s, double h, double lower, double upper) {
  double v_hi = (upper - s) / h;
  double v_lo = (lower - s) / h;
  return epan_F(v_hi) - epan_F(v_lo);
}

// Kernel intensity (log-lambda) for every (subject, grid point) pair,
// replacing the nested lapply/sapply over lambda_fun()/Ker()/edge_c() in
// PLS_Kernel(). Each process is sorted once and swept with a two-pointer
// window against the (assumed ascending) grid, since the kernel has compact
// support [s-h, s+h].
// [[Rcpp::export]]
NumericMatrix log_lambda_matrix_cpp(List PPP_list, NumericVector T_seq,
                                     double h, double lower, double upper) {
  int n = PPP_list.size();
  int T = T_seq.size();
  NumericMatrix out(n, T);

  // Precompute edge correction per grid point (same for every subject).
  std::vector<double> edge(T);
  for (int t = 0; t < T; ++t) edge[t] = edge_c_cpp(T_seq[t], h, lower, upper);

  for (int i = 0; i < n; ++i) {
    NumericVector pp_orig = PPP_list[i];
    std::vector<double> pp(pp_orig.begin(), pp_orig.end());
    std::sort(pp.begin(), pp.end());
    int m = pp.size();

    int lo = 0, hi = 0; // window [lo, hi) of pp within (s-h, s+h)
    for (int t = 0; t < T; ++t) {
      double s = T_seq[t];
      double win_lo = s - h;
      double win_hi = s + h;
      while (lo < m && pp[lo] <= win_lo) ++lo;
      if (hi < lo) hi = lo;
      while (hi < m && pp[hi] < win_hi) ++hi;

      double acc = 0.0;
      for (int k = lo; k < hi; ++k) {
        double v = (s - pp[k]) / h;
        acc += (0.75 / h) * (1.0 - v * v);
      }
      double lambda = acc / edge[t];
      out(i, t) = std::log(lambda);
    }
  }
  return out;
}

// ============================================================================
// Weighted inner product / Gram-Schmidt used to build the PLS basis
// (replaces inprod_mat() and the fully unrolled 10-vector GS()).
// ============================================================================

// [[Rcpp::export]]
double inprod_mat_cpp(const arma::vec& u1, const arma::vec& u2,
                       const arma::mat& K_hat, double lower, double upper) {
  double del = (upper - lower) / u1.n_elem;
  arma::vec tmp = K_hat * u2;
  double value = arma::dot(u1, tmp);
  return del * del * value;
}

// Gram-Schmidt orthonormalization of the columns of U with respect to the
// inner product <u,v> = del^2 * u' K_hat v, generalized to any number of
// columns (the original GS() hardcoded exactly 10).
// [[Rcpp::export]]
arma::mat gram_schmidt_cpp(const arma::mat& U, const arma::mat& K_hat,
                            double lower, double upper) {
  int T = U.n_rows;
  int k = U.n_cols;
  arma::mat V(T, k, arma::fill::zeros);

  for (int j = 0; j < k; ++j) {
    arma::vec v = U.col(j);
    for (int i = 0; i < j; ++i) {
      double c = inprod_mat_cpp(v, V.col(i), K_hat, lower, upper);
      v -= c * V.col(i);
    }
    double nrm = std::sqrt(inprod_mat_cpp(v, v, K_hat, lower, upper));
    V.col(j) = v / nrm;
  }
  return V;
}

// ============================================================================
// Bin counting for the GLM step in P3LS() (replaces the apply() loop over
// bins for every subject).
// ============================================================================

// bins: nb x 2 matrix of (start, end], assumed sorted ascending and
// non-overlapping (as constructed in P3LS()). Returns raw counts (n x nb);
// division by bin length is left to the R caller.
// [[Rcpp::export]]
NumericMatrix bin_counts_matrix_cpp(List PPP_list, NumericMatrix bins) {
  int n = PPP_list.size();
  int nb = bins.nrow();
  NumericMatrix out(n, nb);

  for (int i = 0; i < n; ++i) {
    NumericVector pp_orig = PPP_list[i];
    std::vector<double> pp(pp_orig.begin(), pp_orig.end());
    std::sort(pp.begin(), pp.end());
    int m = pp.size();

    int lo = 0, hi = 0;
    for (int b = 0; b < nb; ++b) {
      double start = bins(b, 0);
      double end = bins(b, 1);
      // advance lo to first index with pp > start
      while (lo < m && pp[lo] <= start) ++lo;
      if (hi < lo) hi = lo;
      while (hi < m && pp[hi] <= end) ++hi;
      out(i, b) = hi - lo;
    }
  }
  return out;
}

// ============================================================================
// Covariance estimator core (replaces the spatstat::crosspairs-based inner
// loop of Cov_estimator() for kernel = "epanechnikov"). Grid is assumed
// sorted ascending (as produced by seq(lbd,ubd,l=ngrid)); the Epanechnikov
// kernel has compact support so only a contiguous window of grid indices
// around each point contributes.
//
// Returns the raw accumulators; the R wrapper combines them with the edge
// correction the same way the original code did:
//   A2   = A_sum + Gsum                    (= sum_i outer(tmp1_i, tmp1_i))
//   Cpld = outer(tmp1_sum, tmp1_sum)        (pooled cross term)
//   A    = A_sum / (n * edge)
//   C    = (Cpld - A2) / (n*(n-1) * edge)
//   R_X  = log(A / C)
// ============================================================================

// [[Rcpp::export]]
List cov_estimator_core_cpp(List PROCESS, NumericVector grids, double bwd) {
  int n = PROCESS.size();
  int ngrid = grids.size();

  arma::mat A_sum(ngrid, ngrid, arma::fill::zeros);
  arma::mat Gsum(ngrid, ngrid, arma::fill::zeros);
  arma::vec tmp1_sum(ngrid, arma::fill::zeros);

  std::vector<double> grid_vec(grids.begin(), grids.end());

  for (int i = 0; i < n; ++i) {
    SEXP proc_sexp = PROCESS[i];
    if (proc_sexp == R_NilValue) continue;
    NumericVector process = proc_sexp;
    int m = process.size();

    arma::vec tmp1(ngrid, arma::fill::zeros);
    arma::mat Gram(ngrid, ngrid, arma::fill::zeros);

    for (int p = 0; p < m; ++p) {
      double x = process[p];
      // window of grid indices with |x - grid[j]| < bwd
      auto lo_it = std::lower_bound(grid_vec.begin(), grid_vec.end(), x - bwd);
      auto hi_it = std::lower_bound(grid_vec.begin(), grid_vec.end(), x + bwd);
      int lo = static_cast<int>(lo_it - grid_vec.begin());
      int hi = static_cast<int>(hi_it - grid_vec.begin());
      if (lo < 0) lo = 0;
      if (hi > ngrid) hi = ngrid;
      if (lo >= hi) continue;

      int wlen = hi - lo;
      std::vector<double> kh(wlen);
      for (int j = 0; j < wlen; ++j) {
        double d = x - grid_vec[lo + j];
        double v = d / bwd; // |v| < 1 within the window by construction
        kh[j] = (0.75 / bwd) * (1.0 - v * v);
        tmp1[lo + j] += kh[j];
      }
      for (int j = 0; j < wlen; ++j) {
        for (int j2 = 0; j2 < wlen; ++j2) {
          Gram(lo + j, lo + j2) += kh[j] * kh[j2];
        }
      }
    }

    A_sum += (tmp1 * tmp1.t()) - Gram;
    Gsum += Gram;
    tmp1_sum += tmp1;
  }

  return List::create(
    Named("A_sum") = A_sum,
    Named("Gsum") = Gsum,
    Named("tmp1_sum") = tmp1_sum
  );
}
