#' Expectation of a function of the follow-up period
#'
#' Integrates a function of the realised follow-up period over the distribution
#' assumed for that period, following Tang (2015, 2018) and Homma and Daimon
#' (2021). In every case \code{dropout} is the proportion of patients who leave
#' before the end of the period.
#'
#' Under \code{"complete"} every patient is followed for the whole period
#' \code{window}. Under \code{"uniform"} a proportion \code{dropout} of patients
#' leave at a time that is uniform on the interval from zero to \code{window}
#' and the remainder are followed to the end, so the period is a two point
#' mixture. Under \code{"exponential"} the period is an exponential variable
#' truncated at \code{window}, with rate \code{-log(1 - dropout) / window}. Both
#' random distributions reduce to complete follow-up when \code{dropout} is
#' zero.
#'
#' @param f Function of a single follow-up period, evaluated on a vector.
#' @param window Prespecified follow-up period.
#' @param dropout Proportion of patients who leave before the end of the period.
#'   A value of one is admissible for the uniform distribution, where it gives
#'   the period used by Tang (2015) when every patient leaves early, but not for
#'   the exponential distribution.
#' @param shape Distribution of the follow-up period, one of
#'   \code{"complete"}, \code{"uniform"} and \code{"exponential"}.
#' @return A numeric value holding the expectation.
#' @importFrom stats integrate
#' @export
#' @examples
#' follow_up_expectation(function(t) t, window = 2, shape = "complete")
#' follow_up_expectation(function(t) t, window = 2, dropout = 0.4,
#'                       shape = "uniform")
#' follow_up_expectation(function(t) t, window = 2, dropout = 0.4,
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
  if (!is.numeric(dropout) || length(dropout) != 1L || is.na(dropout) ||
      dropout < 0 || dropout > 1) {
    stop("'dropout' must be a single numeric value between 0 and 1.")
  }
  if (dropout == 0) {
    return(f(window))
  }
  if (shape == "uniform") {
    early <- integrate(function(t) f(t) / window, lower = 0, upper = window,
                       rel.tol = 1e-10)$value
    return((1 - dropout) * f(window) + dropout * early)
  }
  if (dropout >= 1) {
    stop(paste("'dropout' must be strictly below one when the follow-up",
               "period is exponential."))
  }
  lambda <- -log(1 - dropout) / window
  integrate(function(t) f(t) * lambda * exp(-lambda * t), lower = 0,
            upper = window, rel.tol = 1e-10)$value +
    f(window) * (1 - dropout)
}
