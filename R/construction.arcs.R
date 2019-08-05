construction.arcs <-
function(
  pheromones, capacity, #made in stuart.mmas
  use.order,
  alpha, beta, heuristics
) { #begin function

  #initialize a choice matrix
  solution <- lapply(pheromones,function(x) x<0)

  #initialize chosen vectors
  selected <- lapply(solution,function(x) NA)
  
  for (i in 1:length(capacity)) { #for each factor
    
    #create a pool of possible choices
    pool <- 1:ncol(solution[[i]])

    #randomly select starting position from pool
    selected[[i]] <- sample(pool,1)
    #update pool to exclude starting location
    pool <- pool[!is.element(pool,selected[[i]])]
    #filter out items chosen in other facets
    pool <- pool[!colnames(solution[[i]])[pool]%in%unlist(sapply(solution,function(x) colnames(x)[(colSums(x)+rowSums(x))>0]))]
      

    #filter the special case of one item
    if (capacity[[i]]>1) {
      for (k in 2:capacity[[i]]) {  #for each item
        #compute selection probabilities
        tmp.phe <- pheromones[[i]][selected[[i]][k-1],pool]^alpha
        tmp.heu <- heuristics[[i]][selected[[i]][k-1],pool]^beta
        probs <-  tmp.phe * tmp.heu / sum(tmp.phe*tmp.heu)
        
        #select item (complicated due to sample()-convenience feature)
        if (length(pool)==1) {
          selected[[i]][k] <- pool
        }
        else {
          selected[[i]][k] <- sample(pool,1,FALSE,probs)
        }
        
        #update solution
        solution[[i]][selected[[i]][k-1],selected[[i]][k]] <- TRUE 
        #update pool to exclude choice
        pool <- pool[!is.element(pool,selected[[i]])]
      }        
    }
  }
  
  if (!use.order) {
    selected <- lapply(selected,sort)
  }
  
  return(list(selected=selected,solution=solution))

}
