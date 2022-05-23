#' @export

print.stuartEmpiricalObjective <- function(x, ...) {
  cat('Empirical STUART objective function with:\n\n')
  cat(x$string)
  cat('\n\nUse ...$func() to apply function to data.')
}

#' @export

print.stuartFixedObjective <- function(x, ...) {
  cat('Fixed STUART objective function with:\n\n')
  cat(x$string)
  cat('\n\nUse ...$func() to apply function to data.')
}

#' @export

print.stuartManualObjective <- function(x, ...) {
  cat('Manual STUART objective function with:\n\n')
  cat(x$string)
  cat('\n\nUse ...$func() to apply function to data.')
}