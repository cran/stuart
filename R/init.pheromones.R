init.pheromones <-
function(
  short.factor.structure, localization, alpha
) { #begin function

  pheromones <- list(NA)

  #initialize when depositing on arcs
  if (localization=='arcs') {
    #initialize all pheromones to 1e+100, diagonals to 0
    for (i in 1:length(short.factor.structure)) {
      pheromones[[i]] <- list(NA)
      pheromones[[i]] <- matrix((1e+100)^(1/alpha),length(short.factor.structure[[i]]),length(short.factor.structure[[i]]))
      diag(pheromones[[i]]) <- 0
      dimnames(pheromones[[i]]) <- list(short.factor.structure[[i]],short.factor.structure[[i]])
    }
  }

  #initialize when depositing on nodes
  if (localization=='nodes') {
    #initialize all pheromones to 1e+100
    for (i in 1:length(short.factor.structure)) {
      pheromones[[i]] <- list(NA)
      pheromones[[i]] <- matrix((1e+100)^(1/alpha),1,length(short.factor.structure[[i]]))
      dimnames(pheromones[[i]])[[2]] <- short.factor.structure[[i]]
    }
  }

  names(pheromones) <- names(short.factor.structure)

  return(pheromones)

}
