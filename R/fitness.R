fitness <-
function(objective = NULL, solution.fit
) { #begin function

  criteria <- names(formals(objective))

  output <- list()
  
  if (!all(criteria%in%names(solution.fit))) {
    output[[1]] <- 0
    for (i in 1:length(criteria)) {
      output[[i+1]] <- NA
    }
    names(output) <- c('pheromone',unlist(criteria))
  }
  
  else {
    pheromone <- do.call(objective,solution.fit[names(formals(objective))])
    
    output[[1]] <- pheromone
    for (i in 1:length(criteria)) {
      output[[i+1]] <- solution.fit[[unlist(criteria[i])]]
    }
    names(output) <- c('pheromone',unlist(criteria))
  }
  
  if (length(output$pheromone)!=1) {
    stop('The objective function you provided does not return a single value.', call. = FALSE)
  }
  
  # remove matrices from output
  output$lvcor <- NULL
  output$lambda <- NULL
  output$theta <- NULL
  output$psi <- NULL
  output$alpha <- NULL
  output$beta <- NULL
  
  return(output)

}
