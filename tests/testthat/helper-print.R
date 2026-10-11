quiet_print <- function(x, ...) {
  lines <- capture.output(returned <- withVisible(print(x, ...)))
  list(lines = lines, visible = returned$visible, value = returned$value)
}

base <- list(rate1 = 2.6, rate2 = 1.3, dispersion1 = 0.8, dispersion2 = 0.8,
             shared = 0.4, n1 = 60, n2 = 90, n3 = 300, window1 = 0.5,
             window2 = 1)
