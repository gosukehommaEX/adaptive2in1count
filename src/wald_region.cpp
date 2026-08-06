#include <Rcpp.h>
#include <cmath>
using namespace Rcpp;

//' Rejection region of the Wald test for two negative binomial rates
//'
//' Evaluates the one-sided Wald test on the log rate ratio for every pair of
//' group totals, treating the dispersion as known. A continuity correction of
//' one half is applied to the totals so that the statistic remains finite when
//' a group total is zero. The region is used as a comparator for the exact
//' conditional test.
//'
//' @param size1 Number of patients in the first group.
//' @param size2 Number of patients in the second group.
//' @param dispersion Common dispersion of the negative binomial distribution.
//' @param smax Largest count retained on the support of each group total.
//' @param z Upper standard normal quantile at the one-sided level.
//' @return A numeric matrix of zeros and ones with \code{smax + 1} rows and
//'   \code{smax + 1} columns, indexed by the two group totals.
//' @keywords internal
//' @noRd
// [[Rcpp::export]]
NumericMatrix wald_region_cpp(double size1, double size2, double dispersion, int smax, double z) {
  int S = smax + 1;
  NumericMatrix out(S, S);
  double v1 = 1.0 / (size1 * dispersion);
  double v2 = 1.0 / (size2 * dispersion);
  for (int s1 = 0; s1 < S; s1++) {
    double a1 = (double) s1 + 0.5;
    for (int s2 = 0; s2 < S; s2++) {
      double a2 = (double) s2 + 0.5;
      double lrr = std::log(a1 / size1) - std::log(a2 / size2);
      double se = std::sqrt(1.0 / a1 + v1 + 1.0 / a2 + v2);
      if (lrr / se <= -z) out(s1, s2) = 1.0;
    }
  }
  return out;
}
