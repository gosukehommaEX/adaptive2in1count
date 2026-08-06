#include <Rcpp.h>
using namespace Rcpp;

//' Split of the interim count into its shared and endpoint specific parts
//'
//' The interim count is the sum of a contribution driven by the frailty that
//' the two endpoints share and a contribution specific to the Phase 2 endpoint.
//' Row \code{a} of the returned matrix holds the conditional distribution of the
//' shared contribution given that the interim count equals \code{a}. Rows whose
//' total probability underflows to zero are replaced by a point mass at zero so
//' that every row remains a probability distribution.
//'
//' @param f0 Numeric vector holding the probability mass function of the shared
//'   contribution.
//' @param f1 Numeric vector holding the probability mass function of the
//'   endpoint specific contribution, of the same length as \code{f0}.
//' @return A square numeric matrix whose rows sum to one.
//' @keywords internal
//' @noRd
// [[Rcpp::export]]
NumericMatrix posterior_weight_cpp(NumericVector f0, NumericVector f1) {
  int n = f0.size();
  NumericMatrix out(n, n);
  for (int a = 0; a < n; a++) {
    double total = 0.0;
    for (int x = 0; x <= a; x++) {
      double v = f0[x] * f1[a - x];
      out(a, x) = v;
      total += v;
    }
    if (total > 0.0) {
      for (int x = 0; x <= a; x++) out(a, x) /= total;
    } else {
      out(a, 0) = 1.0;
    }
  }
  return out;
}
