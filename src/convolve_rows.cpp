#include <Rcpp.h>
#include <algorithm>
using namespace Rcpp;

//' Convolve every row of a matrix with a common kernel
//'
//' Each row of \code{x} holds the conditional distribution of the Phase 3 count
//' contributed by the patients who were already followed at the interim
//' analysis. Convolving with \code{kernel}, the distribution of the count
//' contributed by the remaining patients, gives the conditional distribution of
//' the Phase 3 total. The result is truncated at \code{tmax}.
//'
//' @param x Numeric matrix whose rows are probability mass functions.
//' @param kernel Numeric vector holding the probability mass function to
//'   convolve with.
//' @param tmax Largest count retained on the support of the result.
//' @return A numeric matrix with \code{nrow(x)} rows and \code{tmax + 1}
//'   columns.
//' @keywords internal
//' @noRd
// [[Rcpp::export]]
NumericMatrix convolve_rows_cpp(NumericMatrix x, NumericVector kernel, int tmax) {
  int nr = x.nrow();
  int nc = x.ncol();
  int nk = kernel.size();
  NumericMatrix out(nr, tmax + 1);
  for (int i = 0; i < nr; i++) {
    for (int u = 0; u < nc; u++) {
      if (u > tmax) break;
      double xv = x(i, u);
      if (xv == 0.0) continue;
      int vmax = std::min(nk - 1, tmax - u);
      for (int v = 0; v <= vmax; v++) {
        out(i, u + v) += xv * kernel[v];
      }
    }
  }
  return out;
}
