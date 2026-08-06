#include <Rcpp.h>
using namespace Rcpp;

//' Negative binomial probability mass function on a truncated support
//'
//' Evaluates the negative binomial probability mass function with a given
//' dispersion (shape) and mean on the integers 0 to \code{ymax}. The
//' probabilities are accumulated on the log scale so that large shapes do not
//' overflow. A non-positive shape or mean returns a point mass at zero, which
//' is the degenerate limit used by the transition matrix when the shared
//' frailty component or the number of additional patients is zero.
//'
//' @param shape Dispersion (shape) of the negative binomial distribution.
//' @param mean Mean of the negative binomial distribution.
//' @param ymax Largest count retained on the support.
//' @return A numeric vector of length \code{ymax + 1}.
//' @keywords internal
//' @noRd
// [[Rcpp::export]]
NumericVector nb_pmf_cpp(double shape, double mean, int ymax) {
  NumericVector out(ymax + 1);
  if (shape <= 0.0 || mean <= 0.0) {
    out[0] = 1.0;
    return out;
  }
  double p = shape / (shape + mean);
  double lp = std::log(p);
  double lq = std::log1p(-p);
  double lg = 0.0;
  double lf = 0.0;
  for (int y = 0; y <= ymax; y++) {
    if (y > 0) {
      lg += std::log(shape + (double)(y - 1));
      lf += std::log((double) y);
    }
    out[y] = std::exp(lg - lf + shape * lp + (double) y * lq);
  }
  return out;
}

//' Negative binomial probability mass functions for a family of shapes
//'
//' Evaluates one negative binomial probability mass function per element of
//' \code{shapes} and \code{means} on a common support. This is used to build
//' the conditional distribution of the Phase 3 count given the interim count,
//' which is a negative binomial whose shape increases with the interim count.
//'
//' @param shapes Numeric vector of dispersions (shapes).
//' @param means Numeric vector of means, of the same length as \code{shapes}.
//' @param ymax Largest count retained on the support.
//' @return A numeric matrix with one row per element of \code{shapes} and
//'   \code{ymax + 1} columns.
//' @keywords internal
//' @noRd
// [[Rcpp::export]]
NumericMatrix nb_pmf_matrix_cpp(NumericVector shapes, NumericVector means, int ymax) {
  int n = shapes.size();
  NumericMatrix out(n, ymax + 1);
  for (int i = 0; i < n; i++) {
    double shape = shapes[i];
    double mean = means[i];
    if (shape <= 0.0 || mean <= 0.0) {
      out(i, 0) = 1.0;
      continue;
    }
    double p = shape / (shape + mean);
    double lp = std::log(p);
    double lq = std::log1p(-p);
    double lg = 0.0;
    double lf = 0.0;
    for (int y = 0; y <= ymax; y++) {
      if (y > 0) {
        lg += std::log(shape + (double)(y - 1));
        lf += std::log((double) y);
      }
      out(i, y) = std::exp(lg - lf + shape * lp + (double) y * lq);
    }
  }
  return out;
}
