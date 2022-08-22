fitness <-
function(objective = NULL, solution.fit, software, technicals = NULL
) { #begin function

  criteria <- names(formals(objective$func))
  
  output <- list()
  
  if (!all(criteria%in%names(solution.fit))) {
    output[[1]] <- 0
    for (i in 1:length(criteria)) {
      output[[i+1]] <- NA
    }
    names(output) <- c('pheromone',unlist(criteria))
  }
  
  else {
    pheromone <- do.call(objective$func,solution.fit[names(formals(objective$func))])
    
    output[[1]] <- pheromone
    for (i in 1:length(criteria)) {
      output[[i+1]] <- solution.fit[[unlist(criteria[i])]]
    }
    names(output) <- c('pheromone',unlist(criteria))
  }
  
  if (length(output$pheromone)!=1) {
    stop('The objective function you provided does not return a single value.', call. = FALSE)
  }
  
  if (!is.null(technicals)) {
    output$technicals <- solution.fit$technicals[technicals]
  }
  
  # remove matrices from output
  # output$lvcor <- NULL
  # output$lambda <- NULL
  # output$theta <- NULL
  # output$psi <- NULL
  # output$alpha <- NULL
  # output$beta <- NULL
  
  return(output)

}
