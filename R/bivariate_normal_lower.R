#' Lower orthant probability of a standard bivariate normal distribution
#'
#' Evaluates the probability that both coordinates of a standard bivariate
#' normal vector with correlation \code{rho} fall at or below the given
#' quantiles. The probability is obtained by conditioning on the first
#' coordinate and integrating the resulting univariate normal distribution
#' function, so no additional package is required.
#'
#' @param a Upper limit for the first coordinate.
#' @param b Upper limit for the second coordinate.
#' @param rho Correlation between the two coordinates.
#' @return A numeric value holding the probability.
#' @importFrom stats dnorm integrate pnorm
#' @export
#' @examples
#' bivariate_normal_lower(0, 0, rho = 0)
#' bivariate_normal_lower(0, 0, rho = 0.5)
#' bivariate_normal_lower(1.96, -1.96, rho = 0.8)
bivariate_normal_lower <- function(a, b, rho) {
  if (!is.numeric(rho) || length(rho) != 1L || is.na(rho) || abs(rho) > 1) {
    stop("'rho' must be a single numeric value between -1 and 1.")
  }
  if (is.na(a) || is.na(b)) {
    return(NA_real_)
  }
  if (a == -Inf || b == -Inf) {
    return(0)
  }
  if (a == Inf) {
    return(pnorm(b))
  }
  if (b == Inf) {
    return(pnorm(a))
  }
  if (abs(rho) < 1e-12) {
    return(pnorm(a) * pnorm(b))
  }
  if (rho >= 1 - 1e-12) {
    return(pnorm(min(a, b)))
  }
  if (rho <= -1 + 1e-12) {
    return(max(0, pnorm(a) + pnorm(b) - 1))
  }
  # The integrand is negligible outside a few standard deviations, and passing
  # an infinite lower limit together with a distant upper limit makes the
  # transformation used by integrate() miss the region that carries the mass.
  # Both limits are therefore clipped where the standard normal density is
  # below the double precision floor.
  bound <- 8.5
  if (a >= bound) {
    return(pnorm(b))
  }
  if (a <= -bound) {
    return(0)
  }
  root <- sqrt(1 - rho^2)
  value <- integrate(function(u) dnorm(u) * pnorm((b - rho * u) / root),
                     lower = -bound, upper = a, rel.tol = 1e-11)$value
  min(max(value, 0), 1)
}
