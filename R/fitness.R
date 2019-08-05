fitness <-
function(objective = NULL, solution.fit, software
) { #begin function

  criteria <- names(formals(objective))
  
  # translate basic fit components, if necessary
  if (software == 'Mplus') {
    name <- c('rmsea','srmr','wrmr','cfi','tli','chisq','df','pvalue','aic','bic','abic')
    locator <- c('RMSEA_Estimate', 'SRMR', 'WRMR', 'CFI', 'TLI', 'ChiSqM_Value', 'ChiSqM_DF', 'ChiSqM_PValue', 'AIC', 'BIC', 'aBIC')
    
    locator <- locator[which(name%in%criteria)]
    name <- name[which(name%in%criteria)]
    
    for (i in seq_along(locator)) {
      names(solution.fit)[which(names(solution.fit)==locator[i])] <- name[i]
    }
  }

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
