stuart.gene <-
  function(
    short.factor.structure, short, long.equal,    #made on toplevel
    capacity,
    data, factor.structure, auxi, use.order,                       #simple prerequisites
    item.invariance,
    repeated.measures, long.invariance,                            #longitudinal relations
    mtmm, mtmm.invariance,                                         #mtmm relations
    grouping, group.invariance,                                    #grouping relations
    
    software, cores,                                               #Software to be used
    
    objective=NULL, ignore.errors=FALSE,                           #objective function
    
    generations = 128, individuals = 64,                            #settings of the algorithm
    elitism = 1/individuals, reproduction = .5, mutation = .1,
    mating.index = 1, mating.size = .25, 
    mating.criterion = 'fitness',
    tolerance = .0001,
    
    suppress.model=FALSE, analysis.options=NULL,                   #Additional modeling
    seed,
    filename,
    
    ...                                                            #All the other stuff
  ) { #function begin
   
    #set random seed, if provided
    if (!is.null(seed)) {
      old.seed <- .Random.seed
      old.kind <- RNGkind()[1]
      RNGkind("L'Ecuyer-CMRG")
      set.seed(seed)
    }
    
    #initialize fitness results of best solutions
    phe.ib <- 0
    phe.gb <- 0

    #counting
    generation <- 1
        
    ### Loops ###
    log <- list()
    conv <- NULL
    
    #creating user feedback
    message('Running STUART with Genetic Algorithm.\n')
    progress <- utils::txtProgressBar(0,max(c(generations,1)),style=3)
    
    # generate initial population
    full <- FALSE
    n <- individuals
    combinations <- do.call('generate.combinations',mget(names(formals(generate.combinations))))
    
    duplicate <- combinations$duplicate
    filter <- combinations$filter[!duplicated(duplicate), , drop = FALSE]
    combi <- combinations$combi
    tried <- matrix(NA, ncol = sum(unlist(capacity)))[-1,]
    
    if (software=='Mplus') {
      #file location
      if (is.null(filename)) filename <- paste0(tempdir(), '/stuart')
      
      #writing the data file
      utils::write.table(data,paste(filename,'_data.dat',sep=''),
        col.names=FALSE,row.names=FALSE,na='-9999',
        sep='\t',dec='.')
    }
    
    repeat { #over generations
      
      output.model <- FALSE
      svalues <- FALSE
      bf.args <- mget(names(formals(bf.cycle))[-1])
      combi_mat <- do.call(cbind, combi)
      
      if (nrow(filter) > 0) {
        #parallel processing for R-internal estimations
        if (software=='lavaan') {
          if (cores>1) {
            #set up parallel processing on windows
            if (grepl('Windows',Sys.info()[1],ignore.case=TRUE)) {
              cl <- parallel::makeCluster(cores)
              
              bf.results <- parallel::parLapply(cl,1:nrow(filter),function(run) {
                do.call('bf.cycle',c(run,bf.args))
              })
              parallel::stopCluster(cl)
            }
            
            #run ants in parallel on unixies
            else {
              bf.results <- parallel::mclapply(1:nrow(filter), function(run) {
                do.call('bf.cycle',c(run,bf.args))},
                mc.cores=cores
              )
            }
          } else {
            bf.results <- lapply(1:nrow(filter),function(run) {
              do.call('bf.cycle',c(run,bf.args))})
          }
        }
        
        #serial processing if Mplus is used (Mplus-internal parallelization is used)
        if (software=='Mplus') {
          bf.args$filename <- filename
          bf.args$cores <- cores
          bf.results <- lapply(1:nrow(filter), function(run) {     
            do.call('bf.cycle',c(run,bf.args))})
        }
      }
      
      #fill in results for duplicates
      tmp <- vector('list', individuals)
      tmp[filter[,1]] <- bf.results
      if (generation == 1) {
        bf.results <- tmp[duplicate]
      } else {
        tmp[sapply(tmp,is.null)] <- log[stats::na.omit(duplicate)]
        bf.results <- tmp
      }
      
      # components of genetic algorithm
      # parent selection: fitness proportionate selection
      pheromones <- sapply(bf.results, function(x) x$solution.phe$pheromone)
      if (sum(pheromones > 0) == 0) {
        stop(paste0('There were no viable solutions in generation ', generation,'. This may indicate estimation problems in the CFA.'), call. = FALSE)
      }
      if (sum(pheromones > 0) < round(individuals * reproduction)) {
        warning(paste0('There were not enough viable individuals in generation ', generation, '. Some non-viables were randomly selected.\n'), call. = FALSE)
        parents <- (1:individuals)[pheromones > 0]
        parents <- c(parents, sample((1:individuals)[-parents], round(individuals * reproduction)-length(parents)))
      }
      else {
        parents <- sample(1:individuals, round(individuals * reproduction), FALSE, pheromones / sum(pheromones))
      }
      
      # random mating
      if (is.null(mating.index)) {
        mating <- matrix(sample(rep_len(parents, individuals*2)), ncol = 2)
      } 
      # fitness and similarity mating
      else {
        parents <- rep_len(parents, individuals)
        mating <- matrix(NA, ncol = 2, nrow = individuals)
        for (i in seq_along(parents)) {
          partners <- sample(parents[-i], max(individuals*reproduction*mating.size, 1))
          if (mating.criterion == 'fitness') {
            mating[i,] <- c(parents[i], 
              partners[which(pheromones[partners] == 
                  stats::quantile(pheromones[partners], mating.index, type = 1))[1]])
          }
          if (mating.criterion == 'similarity') {
            tmp <- combi_mat[partners,]
            similar <- apply(tmp, 1, function(x) length(intersect(combi_mat[parents[i],], x))/length(union(combi_mat[parents[i],], x)))
            mating[i,] <- c(parents[i], 
              partners[which(similar == stats::quantile(similar, mating.index, type = 1))[1]])
          }
        }
      }
      
      # elitism
      nextgen <- lapply(combi, function(x) x[which(rank(pheromones, ties.method = 'random')>round(individuals*(1-elitism))), , drop = FALSE])
      if (nrow(nextgen[[1]]) > 0) {
        nextgen <- lapply(nextgen, function(x) x[!duplicated(do.call(cbind, nextgen)), , drop = FALSE])
      }

      for (i in 1:nrow(mating)) {
        dad <- lapply(combi, function(x) x[mating[i,1],])
        mom <- lapply(combi, function(x) x[mating[i,2],])

        for (i in 1:length(short.factor.structure)) {
          dad_solution <- seq_along(short.factor.structure[[i]])%in%dad[[i]]
          mom_solution <- seq_along(short.factor.structure[[i]])%in%mom[[i]]
          dif_solution <- dad_solution - mom_solution
          tmp <- which(cumsum(dif_solution) == 0)
          crossover <- ifelse(length(tmp)>1, sample(tmp, 1), tmp)
          if (is.na(crossover) | crossover == length(dad_solution)) crossover <- 0
          
          kid_solution <- c(dad_solution[0:crossover],mom_solution[(crossover+1):length(mom_solution)])
          if (sample(c(TRUE,FALSE), 1, FALSE, c(mutation, 1-mutation))) {
            tmp <- c(sample(which(kid_solution), 1),sample(which(!kid_solution), 1))
            kid_solution[tmp] <- kid_solution[rev(tmp)]
          }
          nextgen[[i]] <- rbind(nextgen[[i]], which(kid_solution))
        }
      }
      
      nextgen <- lapply(nextgen, function(x) x[1:individuals, ])
      
      #iteration.best memory
      individual.ib <- which.max(sapply(bf.results, function(x) return(x$solution.phe$pheromone)))
      solution.ib <- list()
      for (i in seq_along(short.factor.structure)) {
        solution.ib[[i]] <- seq_along(short.factor.structure[[i]])%in%bf.results[[individual.ib]]$selected[[i]]
      }
      phe.ib <- bf.results[[individual.ib]]$solution.phe$pheromone
      selected.ib <- bf.results[[individual.ib]]$selected
      
      #global.best memory
      if (phe.ib > phe.gb | generation == 1) {
        solution.gb <- solution.ib
        phe.gb <- phe.ib
        selected.gb <- selected.ib
      }
      
      #create matrix of combinations already evaluated
      tried <- rbind(tried, combi_mat)
      
      #check for duplication in nextgen
      duplicate <- match(data.frame(t(do.call(cbind, nextgen))), data.frame(t(tried)))
      filter <- data.frame(matrix(which(is.na(duplicate)),
        nrow=sum(is.na(duplicate)),
        ncol=length(short.factor.structure)))
      
      #create log
      log <- c(log, bf.results)
      utils::setTxtProgressBar(progress,generation)
      
      # check for convergence
      conv <- c(conv, phe.gb)
      if (generation > max(min(c(.1*generations, 10)),1) & stats::var(conv/conv[1]) <= tolerance) {
        end.reason <- 'Algorithm converged.'
        break
      }
      if (generation >= generations) {
        end.reason <- 'Maximum number of generations exceeded.'
        break
      }
      else {
        #go on to next generation
        combi <- nextgen
        generation <- generation + 1
      }
      
    } # end loop
     
    #feedback
    close(progress)
    message(paste('\nSearch ended.',end.reason))      
    
    # reformat log
    log <- cbind(rep(1:generation,each = individuals),rep(1:individuals, generation),t(sapply(log, function(x) array(data=unlist(x$solution.phe)))))
    log <- data.frame(log)
    names(log) <- c('run','ind',names(bf.results[[1]]$solution.phe))
    
    
    #return to previous random seeds
    if (!is.null(seed)) {
      RNGkind(old.kind)
      .Random.seed <<- old.seed
    }
    
    #Generating Output
    for (i in seq_along(solution.gb)) names(solution.gb[[i]]) <- short.factor.structure[[i]]
    results <- mget(grep('.gb',ls(),value=TRUE))
    results$selected.items <- translate.selection(selected.gb,factor.structure,short)
    results$log <- log
    results$pheromones <- pheromones
    results$parameters <- list(generations, individuals, elitism, mutation, mating.index, mating.size, 
      mating.criterion, tolerance, var.gb = stats::var(conv/conv[1]),
      seed, objective, factor.structure)
    names(results$parameters) <- c('generations', 'individuals', 'elitism', 'mutation', 
      'mating.index', 'mating.size', 'mating.criterion', 'tolerance', 'var.gb',
      'seed', 'objective', 'factor.structure')
    return(results)
    
  }