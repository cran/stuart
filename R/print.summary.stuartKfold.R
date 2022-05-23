#' @export

print.summary.stuartKfold <-
  function(x,...) {
    message('Warning: This is a beta-build of stuart. Please report any bugs you encounter.\n')
    cat('SUMMARY OF ANALYSIS:\n\n')
    cat('Number of Folds:', x$k, '\n')
    cat('Analysis Type:', x$Type, '\n')
    cat('Estimation Software:', x$Software, '\n')
    cat('Models Estimated:', x$Models, '\n')
    cat('Time Required:',x$Time,'seconds\n')
    cat('\nCrossvalidation Results with', toupper(x$max.invariance), 'Invariance:\n')
    print(x$Results)
    cat('\nAverage Jaccard Similarity: ')
    cat(paste0(paste0(names(x$Jaccard),': ',paste(x$Jaccard), collapse=', '), '\n'))
    cat(paste0('\nConstructed Subtests: (k = ', which.max(x$Results$pheromone), ')\n'))
    for (i in 1:length(x$Subtests)) {
      cat(paste0(names(x$Subtests)[i],': ',paste(x$Subtests[[i]],collapse=' '),'\n'))
    }
    
  }
