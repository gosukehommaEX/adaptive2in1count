#include <Rcpp.h>
#include <algorithm>
using namespace Rcpp;

//' Distribution of the Phase 2 total given the interim count
//'
//' The Phase 2 total is the interim count plus an independent contribution from
//' the patients enrolled after the interim analysis. Row \code{a} of the
//' returned matrix therefore holds \code{f} shifted by \code{a} positions and
//' truncated at \code{smax}.
//'
//' @param f Numeric vector holding the probability mass function of the
//'   contribution from the patients enrolled after the interim analysis.
//' @param amax Largest count retained on the support of the interim count.
//' @param smax Largest count retained on the support of the Phase 2 total.
//' @return A numeric matrix with \code{amax + 1} rows and \code{smax + 1}
//'   columns.
//' @keywords internal
//' @noRd
// [[Rcpp::export]]
NumericMatrix shift_matrix_cpp(NumericVector f, int amax, int smax) {
  NumericMatrix out(amax + 1, smax + 1);
  int nf = f.size();
  for (int a = 0; a <= amax; a++) {
    if (a > smax) break;
    int bmax = std::min(nf - 1, smax - a);
    for (int b = 0; b <= bmax; b++) out(a, a + b) = f[b];
  }
  return out;
}
