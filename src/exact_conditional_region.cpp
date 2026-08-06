#include <Rcpp.h>
#include <vector>
#include <algorithm>
using namespace Rcpp;

//' Rejection region of the exact conditional test for two negative binomial rates
//'
//' Under the null hypothesis of equal rates, equal exposure per patient and a
//' known common dispersion, the two group totals are negative binomial with a
//' common success probability. Conditioning on their sum removes that
//' probability and leaves a negative hypergeometric distribution, so the
//' resulting test is exact and free of the nuisance rate. The function returns
//' the indicator of the one-sided rejection region in favour of the first group
//' having the lower rate.
//'
//' @param shape1 Dispersion times the number of patients in the first group.
//' @param shape2 Dispersion times the number of patients in the second group.
//' @param smax Largest count retained on the support of each group total.
//' @param alpha One-sided significance level.
//' @return A numeric matrix of zeros and ones with \code{smax + 1} rows and
//'   \code{smax + 1} columns, indexed by the two group totals.
//' @keywords internal
//' @noRd
// [[Rcpp::export]]
NumericMatrix exact_conditional_region_cpp(double shape1, double shape2, int smax, double alpha) {
  int S = smax + 1;
  NumericMatrix out(S, S);
  std::vector<double> lc1(S), lc2(S);
  double lg1 = R::lgammafn(shape1);
  double lg2 = R::lgammafn(shape2);
  for (int u = 0; u < S; u++) {
    lc1[u] = R::lgammafn((double) u + shape1) - R::lgammafn((double) u + 1.0) - lg1;
    lc2[u] = R::lgammafn((double) u + shape2) - R::lgammafn((double) u + 1.0) - lg2;
  }
  int tmax = 2 * smax;
  std::vector<double> lct(tmax + 1);
  double lgs = R::lgammafn(shape1 + shape2);
  for (int t = 0; t <= tmax; t++) {
    lct[t] = R::lgammafn((double) t + shape1 + shape2) - R::lgammafn((double) t + 1.0) - lgs;
  }
  std::vector<double> w(S);
  for (int t = 0; t <= tmax; t++) {
    int lo = std::max(0, t - smax);
    int hi = std::min(smax, t);
    double total = 0.0;
    for (int u = lo; u <= hi; u++) {
      double v = std::exp(lc1[u] + lc2[t - u] - lct[t]);
      w[u - lo] = v;
      total += v;
    }
    if (total <= 0.0) continue;
    double cum = 0.0;
    for (int u = lo; u <= hi; u++) {
      cum += w[u - lo] / total;
      if (cum <= alpha + 1e-12) {
        out(u, t - u) = 1.0;
      } else {
        break;
      }
    }
  }
  return out;
}
