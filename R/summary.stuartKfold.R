#' @export

summary.stuartKfold <- function(object, ...) {
  filt <- sapply(object$full, function(x) length(x) == 0) | sapply(object$crossvalidations, function(x) length(x) == 0)
  Time <- as.numeric(sum(sapply(object$full[!filt], function(x) x$timer[3])))
  Models <- sum(sapply(object$full[!filt], function(x) nrow(x$log)))
  Type <- object$call$type
  k <- object$call$k
  Software <- ifelse(is.null(object$call$software), 'lavaan', object$call$software)
  max.invariance <- ifelse(is.null(object$call$max.invariance), 'strict', object$call$max.invariance)
  
  Results <- do.call('rbind', lapply(object$crossvalidations[!filt], function(x) x[nrow(x), ]))
  rownames(Results) <- paste('k =', 1:k)[!filt]
  
  # Compute Jaccard index
  tmp <- lapply(object$full[!filt], `[[`, 'solution')
  solu <- list()
  for (i in seq_along(tmp[[1]])) {
    solu[[i]] <- do.call('rbind', lapply(tmp, function(x) x[[i]]))
  }
  names(solu) <- names(tmp[[1]])
  Jaccard <- lapply(lapply(lapply(solu, stats::dist, method = 'binary'), mean, na.rm = TRUE), function(x) 1-x)

  
  Subtests <- object$subtests
  
  out <- list(Time, Models, Type, k, Software, max.invariance, Results, Jaccard, Subtests)
  names(out) <- c('Time', 'Models', 'Type', 'k', 'Software', 'max.invariance', 'Results', 'Jaccard', 'Subtests')
  class(out) <- 'summary.stuartKfold'
  return(out)
}