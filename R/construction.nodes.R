construction.nodes <-
function(
  pheromones, capacity, #made in stuart.mmas
  use.order,
  alpha, beta, heuristics  
) { #begin function

  solution <- lapply(pheromones,function(x) x<0)

  #initialize chosen vectors
  selected <- lapply(solution,function(x) NA)

  for (i in 1:length(capacity)) { #for each factor
    #create a pool of possible choices
    pool <- 1:ncol(solution[[i]])
    #filter out items chosen in other facets
    pool <- pool[!colnames(solution[[i]])[pool]%in%unlist(sapply(solution,function(x) colnames(x)[x]))]
    
    #compute selection probabilities
    tmp.phe <- pheromones[[i]][pool]^alpha
    tmp.heu <- heuristics[[i]][pool]^beta
    probs <-  tmp.phe * tmp.heu / sum(tmp.phe*tmp.heu)

    #select items
    selected[[i]] <- sample(pool,capacity[[i]],FALSE,probs)

    solution[[i]][selected[[i]]] <- TRUE
  }

  if (!use.order) {
    selected <- lapply(selected,function(x) sort(x))
  }

  return(list(selected=selected,solution=solution))

}
