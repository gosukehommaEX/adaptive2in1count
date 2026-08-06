#' Expand a per-group sample size to a vector of length two
#'
#' Accepts either a single sample size, which is used for both groups, or a
#' vector of length two giving the investigational and the control group, and
#' returns a validated vector of length two.
#'
#' @param size Numeric vector of length one or two holding the sample size.
#' @param name Name of the argument, used in the error message.
#' @return A numeric vector of length two.
#' @export
#' @examples
#' expand_group_size(50, "n1")
#' expand_group_size(c(60, 40), "n1")
expand_group_size <- function(size, name = "size") {
  if (!is.numeric(size) || anyNA(size) || any(size <= 0) ||
      any(size != floor(size)) || !(length(size) %in% c(1L, 2L))) {
    stop(sprintf("'%s' must be one or two positive whole numbers.", name))
  }
  if (length(size) == 1L) {
    size <- rep(size, 2L)
  }
  as.numeric(size)
}
