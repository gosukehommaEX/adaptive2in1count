#' Expectation of a function of the follow-up period
#'
#' Integrates a function of the realised follow-up period over the distribution
#' assumed for that period. Three distributions are available, following Tang
#' (2015, 2018) and Homma and Daimon (2021). Under \code{"complete"} every
#' patient is followed for the whole period \code{window}. Under
#' \code{"uniform"} the follow-up period is uniform on the interval from zero to
#' \code{window}, and \code{dropout} is not used. Under \code{"exponential"} the
#' follow-up period is an exponential variable truncated at \code{window}, with
#' rate \code{-log(1 - dropout) / window}, so that \code{dropout} is the
#' proportion of patients who leave before the end of the period.
#'
#' @param f Function of a single follow-up period, evaluated on a vector.
#' @param window Prespecified follow-up period.
#' @param dropout Proportion of patients who drop out before the end of the
#'   period, used only by the exponential distribution.
#' @param shape Distribution of the follow-up period, one of
#'   \code{"complete"}, \code{"uniform"} and \code{"exponential"}.
#' @return A numeric value holding the expectation.
#' @importFrom stats integrate
#' @export
#' @examples
#' follow_up_expectation(function(t) t, window = 2, shape = "complete")
#' follow_up_expectation(function(t) t, window = 2, shape = "uniform")
#' follow_up_expectation(function(t) t, window = 2, dropout = 0.3,
#'                       shape = "exponential")
follow_up_expectation <- function(f, window, dropout = 0,
                                  shape = c("complete", "uniform",
                                            "exponential")) {
  shape <- match.arg(shape)
  if (!is.numeric(window) || length(window) != 1L || is.na(window) ||
      window <= 0) {
    stop("'window' must be a single positive numeric value.")
  }
  if (shape == "complete") {
    return(f(window))
  }
  if (shape == "uniform") {
    return(integrate(function(t) f(t) / window, lower = 0, upper = window,
                     rel.tol = 1e-10)$value)
  }
  if (!is.numeric(dropout) || length(dropout) != 1L || is.na(dropout) ||
      dropout <= 0 || dropout >= 1) {
    stop(paste("'dropout' must be a single numeric value strictly between 0",
               "and 1 when the follow-up period is exponential."))
  }
  lambda <- -log(1 - dropout) / window
  integrate(function(t) f(t) * lambda * exp(-lambda * t), lower = 0,
            upper = window, rel.tol = 1e-10)$value +
    f(window) * (1 - dropout)
}
