#' @export

print.summary.stuartOutput <-
function(x,...) {
  message('Warning: This is a beta-build of stuart. Please report any bugs you encounter.\n')
  cat('SUMMARY OF ANALYSIS:\n\n')
  cat('Analysis Type:',x$Type[1],'\n')
  cat('Estimation Software:',x$Software,'\n')
  cat('Models Estimated:',x$Models,'\n')
  cat('Replications of final solution:',x$Replications,'\n')
  if ('end.reason' %in% names(x)) cat(x$end.reason, '\n')
  cat('Time Required:',x$Time,'seconds\n')
  cat('\nOptimization History:\n')
  print(x$Results, row.names = FALSE)
  cat('\nConstructed Subtests:\n')
  paste(paste(names(x$Subtests),': ',paste(x$Subtests,collapse=' '),'\n',collapse=' ',sep=''))
  for (i in 1:length(x$Subtests)) {
    cat(paste0(names(x$Subtests)[i],': ',paste(x$Subtests[[i]],collapse=' '),'\n'))
  }

}
