manualobjective <- function(x, ...) {
  if (!inherits(x, 'function')) {
    stop('The manual objective you provided is not a function.', call. = FALSE)
  }
  arguments <- names(formals(x))
  body <- toString(body(x)[-1])
  add <- vector('character')
  criteria <- vector('character')
  for (i in arguments) {
    filt <- grepl(i, body)
    if (filt) {
      criteria <- c(criteria, i)
    } else {
      add <- c(add, i)
    }
  }
  string <- body
  func <- x
  called <- list(criteria = criteria, add = add)
  out <- list(func = func, string = string, call = called)
  
  class(out) <- 'stuartManualObjective'
  return(out)
}
