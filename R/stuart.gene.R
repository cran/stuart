stuart.gene <-
  function(
    short.factor.structure, short, long.equal, comparisons.equal,
    comparisons.invariance, #made on toplevel
    capacity,
    data, factor.structure, auxi, use.order,                       #simple prerequisites
    item.invariance,
    repeated.measures, long.invariance,                            #longitudinal relations
    mtmm, mtmm.invariance,                                         #mtmm relations
    grouping, group.invariance,                                    #grouping relations
    comparisons,
    software, cores,                                               #Software to be used
    
    objective=NULL, ignore.errors=FALSE, burnin = 5,               #objective function
    
    generations = 256, individuals = 64,                                  #algorithm specs
    selection = 'tournament', selection.pressure = NULL,
    elitism = NULL, reproduction = .5, mutation = .05,
    mating.index = 0, mating.size = .25, 
    mating.criterion = 'similarity',
    immigration = 0,
    convergence.criterion = 'geno.between',
    tolerance = NULL,
    
    reinit.n = 1, reinit.criterion = convergence.criterion,
    reinit.tolerance = NULL, reinit.prop = .75,
    
    schedule = 'run',
    
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

    
    # bind dependent parameters together
    if (is.null(selection.pressure)) {
      selection.pressure <- selection
      selection.pressure[selection == 'tournament'] <- 5
      selection.pressure[selection == 'proportional'] <- 1
      if (is.matrix(selection.pressure)) {
        selection.pressure <- matrix(as.numeric(selection.pressure), ncol = ncol(selection.pressure))
      } else {
        selection.pressure <- as.numeric(selection.pressure)
      }
    }
    
    if (is.null(elitism)) {
      elitism <- individuals
      if (is.matrix(elitism)) {
        elitism[, 2] <- 1/elitism[, 2]
      } else {
        elitism <- 1/elitism
      }
    }
    
    # sanity check for tolerance-convergence.criterion compatibility
    if ((!is.null(tolerance) & !is.list(tolerance) & (length(convergence.criterion) > 1)) |
        (!is.null(reinit.tolerance) & !is.list(reinit.tolerance) & (length(reinit.criterion) > 1))) {
      stop(paste0('When using multiple convergence or reinitialization criteria, the accompanying tolerance must be either NULL or a list.'), .call = FALSE)
    }
    
    # decompress tolerances (for multiple convergence criteria)
    # explicit assignment to avoid check note
    if (is.null(tolerance)) tolerance <- vector('list', length(convergence.criterion))
    if (!is.list(tolerance)) tolerance <- list(tolerance)
    if (is.null(names(tolerance))) names(tolerance) <- convergence.criterion
    tolerance_va <- tolerance[['variance']]
    tolerance_md <- tolerance[['median']]
    tolerance_gw <- tolerance[['geno.within']]
    tolerance_gb <- tolerance[['geno.between']]

    if (is.null(reinit.tolerance)) reinit.tolerance <- vector('list', length(reinit.criterion))
    if (!is.list(reinit.tolerance)) reinit.tolerance <- list(reinit.tolerance)
    if (is.null(names(reinit.tolerance))) names(reinit.tolerance) <- reinit.criterion
    reinit.tolerance_va <- reinit.tolerance[['variance']]
    reinit.tolerance_md <- reinit.tolerance[['median']]
    reinit.tolerance_gw <- reinit.tolerance[['geno.within']]
    reinit.tolerance_gb <- reinit.tolerance[['geno.between']]
    
    # tolerance presets
    tols <- ls(pattern = 'tolerance_')
    pres <- c(.05, .7, .10, .005,
      .01, .8, .05, .0005)
    for (i in seq_along(tols)) {
      if (is.null(get(tols[i]))) assign(tols[i], pres[i])
    }

    #initialize scheduling 
    scheduled <- c('generations', 'individuals', 'selection', 'selection.pressure',
      'elitism', 'reproduction', 'mutation', 'mating.index', 'mating.size', 'mating.criterion',
      'immigration', 'tolerance_va', 'tolerance_md', 'tolerance_gw', 'tolerance_gb', 'reinit.n', 'reinit.criterion', 'reinit.tolerance_va', 'reinit.tolerance_md', 'reinit.tolerance_gw', 'reinit.tolerance_gb', 'reinit.prop')
    
    #global assignment to avoid check note
    generations_cur <- NA 
    individuals_cur <- NA 
    selection_cur <- NA 
    selection.pressure_cur <- NA 
    elitism_cur <- NA 
    reproduction_cur <- NA 
    mutation_cur <- NA 
    mating.index_cur <- NA 
    mating.size_cur <- NA 
    mating.criterion_cur <- NA 
    immigration_cur <- NA 
    tolerance_va_cur <- NA 
    tolerance_md_cur <- NA 
    tolerance_gw_cur <- NA 
    tolerance_gb_cur <- NA 
    reinit.n_cur <- NA 
    reinit.criterion_cur <- NA 
    reinit.tolerance_va_cur <- NA 
    reinit.tolerance_md_cur <- NA 
    reinit.tolerance_gw_cur <- NA 
    reinit.tolerance_gb_cur <- NA 
    reinit.prop_cur <- NA
    
    
    filt <- sapply(mget(scheduled),is.array)
    for (i in 1:length(scheduled[!filt])) {
      assign(paste0(scheduled[!filt][i],'_cur'),mget(scheduled[!filt][i])[[1]])
    }
    if (length(scheduled[filt])>0) {
      scheduled <- scheduled[filt]
      for (i in 1:length(scheduled)) {
        tmp <- mget(scheduled[i])[[1]]
        if (!any(c(0,1)%in%tmp[,1])) {
          stop(paste('The parameter schedule for',scheduled[i],'does not contain a value for the first generation.'),call.=FALSE)
        }
        tmp <- tmp[which.min(tmp[,1]),2]
        assign(paste0(scheduled[i],'_cur'),tmp)
      }
    } else {
      scheduled <- NULL
    }
    
    #counting
    generation <- 1
    run <- 1
    
    # reinitialization indicator
    reinit <- FALSE
    
    # genotypes
    geno <- list()
        
    ### Loops ###
    log <- list()
    qual.ib <- NULL
    cur_reinit.n <- reinit.n_cur
    
    #creating user feedback
    message('Running STUART with Genetic Algorithm.\n')
    progress <- utils::txtProgressBar(0,max(c(generations_cur,1)),style=3)
    
    # generate initial population
    full <- FALSE
    n <- individuals_cur
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
      tmp <- vector('list', individuals_cur)
      tmp[filter[,1]] <- bf.results
      if (run == 1) {
        bf.results <- tmp[duplicate]
      } else {
        redo <- lapply(log[stats::na.omit(duplicate)], function(x) {
          if(all(is.na(x$solution.phe[-1]))) x$solution.phe$pheromone <- 0
          else x$solution.phe$pheromone <- do.call(objective$func, x$solution.phe[-1])
          if(is.na(x$solution.phe$pheromone)) x$solution.phe$pheromone <- 0
          return(x)})
        tmp[sapply(tmp,is.null)] <- redo
        bf.results <- tmp
      }
      
      #parameter schedule
      if (!is.null(scheduled)) {
        for (i in 1:length(scheduled)) {
          tmp <- mget(scheduled[i])[[1]]
          if (schedule=='run') {
            if (any(tmp[,1]==run)) {
              message(paste0('Scheduled value of ',scheduled[i],' updated to ',tmp[which(tmp[,1]==run),2],'.'))
            }
            tmp <- tmp[max(which(tmp[,1]<=run)),2]
          } 
          if (schedule=='generation') {
            if (any(tmp[,1]==generation)) {
              message(paste0('Scheduled value of ',scheduled[i],' updated to ',tmp[which(tmp[,1]==generation),2],'.'))
            }
            tmp <- tmp[max(which(tmp[,1]<=generation)),2]
          }
          assign(paste0(scheduled[i],'_cur'),tmp)
        }
      }
      
      cur_reinit.n <- min(cur_reinit.n, reinit.n_cur)
      
      # components of genetic algorithm
      # parent selection: fitness proportionate selection
      pheromones <- sapply(bf.results, function(x) x$solution.phe$pheromone)
      if (sum(pheromones > 0) == 0) {
        stop(paste0('There were no viable solutions in generation ', generation,'. This may indicate estimation problems in the CFA.'), call. = FALSE)
      }
      if (sum(pheromones > 0) < round(individuals_cur * reproduction_cur) & 
          length(pheromones) >= round(individuals_cur * reproduction_cur)) {
        warning(paste0('There were not enough viable individuals in generation ', generation, '. Some non-viables were randomly selected.\n'), call. = FALSE)
        parents <- (1:length(pheromones))[pheromones > 0]
        parents <- c(parents, sample((1:individuals_cur)[-parents], round(individuals_cur * reproduction_cur)-length(parents)))
      }
      else {
        if (!(selection_cur %in% c('proportional', 'tournament'))) {
          stop(paste0('The selection must be either proportional or tournament. You provided ', selection_cur, '.'), call. = FALSE)
        }
        if (selection_cur == 'proportional') {
          if (length(pheromones) >= round(individuals_cur * reproduction_cur)) {
            parents <- sample(1:length(pheromones), round(individuals_cur * reproduction_cur), FALSE, pheromones^selection.pressure_cur / sum(pheromones^selection.pressure_cur))            
          } else {
            tmp <- floor(round(individuals_cur * reproduction_cur) / length(pheromones))
            parents <- rep(1:length(pheromones), tmp)
            parents <- c(parents, sample(1:length(pheromones), round(individuals_cur * reproduction_cur) - length(parents), FALSE, pheromones^selection.pressure_cur / sum(pheromones^selection.pressure_cur)))
          }
        }
        if (selection_cur == 'tournament') {
          parents <- rep(NA, round(individuals_cur * reproduction_cur))
          pool <- 1:length(pheromones)
          for (i in 1:round(individuals_cur * reproduction_cur)) {
            if (length(pool) == 0) pool <- 1:length(pheromones)
            if (length(pool) < selection.pressure_cur) {
              tmp <- pool
            } else {
              tmp <- sample(pool, selection.pressure_cur)
            }
            parents[i] <- tmp[which.max(pheromones[tmp])]
            pool <- pool[pool!=parents[i]]
          }
        }
      }
      
      # random mating
      if (is.null(mating.index_cur)) {
        mating <- matrix(sample(rep_len(parents, individuals_cur*2)), ncol = 2)
      } 
      # fitness and similarity mating
      else {
        parents <- rep_len(parents, individuals_cur)
        mating <- matrix(NA, ncol = 2, nrow = individuals_cur)
        for (i in seq_along(parents)) {
          partners <- sample(parents[-i], max(individuals_cur*reproduction_cur*mating.size_cur, 1))
          if (mating.criterion_cur == 'fitness') {
            mating[i,] <- c(parents[i], 
              partners[which(pheromones[partners] == 
                  stats::quantile(pheromones[partners], mating.index_cur, type = 1))[1]])
          }
          if (mating.criterion_cur == 'similarity') {
            tmp <- combi_mat[partners,,drop=FALSE]
            similar <- apply(tmp, 1, function(x) length(intersect(combi_mat[parents[i],], x))/length(union(combi_mat[parents[i],], x)))
            mating[i,] <- c(parents[i], 
              partners[which(similar == stats::quantile(similar, mating.index_cur, type = 1))[1]])
          }
        }
      }
      
      # elitism
      nextgen <- lapply(combi, function(x) x[which(rank(pheromones, ties.method = 'random')>round(individuals_cur*(1-elitism_cur))), , drop = FALSE])
      if (nrow(nextgen[[1]]) > 0) {
        nextgen <- lapply(nextgen, function(x) x[!duplicated(do.call(cbind, nextgen)), , drop = FALSE])
      }
      
      # immigration
      if (immigration_cur > 0) {
        n <- round(individuals_cur * immigration_cur)
        combinations <- do.call('generate.combinations',mget(names(formals(generate.combinations))))
        for (i in 1:length(factor.structure)) {
          nextgen[[i]] <- rbind(nextgen[[i]], combinations$combi[[i]])
        }
      }
      
      # mating
      for (i in 1:nrow(mating)) {
        dad <- lapply(combi, function(x) x[mating[i,1],])
        mom <- lapply(combi, function(x) x[mating[i,2],])

        for (i in 1:length(short.factor.structure)) {
          dad_solution <- seq_along(short.factor.structure[[i]])%in%dad[[i]]
          mom_solution <- seq_along(short.factor.structure[[i]])%in%mom[[i]]
          dif_solution <- dad_solution - mom_solution
          tmp <- which(cumsum(dif_solution) == 0)
          crossover <- ifelse(length(tmp)>1, sample(tmp, 1), tmp)
          if (is.na(crossover) | crossover == length(dad_solution)) {
            kid_solution <- get(sample(c('mom_solution', 'dad_solution'), 1))
          } else {
            kid_solution <- c(dad_solution[0:crossover],mom_solution[(crossover+1):length(mom_solution)])
          }

          if (sample(c(TRUE,FALSE), 1, FALSE, c(mutation_cur, 1-mutation_cur))) {
            tmp <- c(sample(which(kid_solution), 1), sample(which(!kid_solution), 1))
            kid_solution[tmp] <- kid_solution[rev(tmp)]
          }
          nextgen[[i]] <- rbind(nextgen[[i]], which(kid_solution))
        }
      }
      
      nextgen <- lapply(nextgen, function(x) x[1:individuals_cur, ])
      
      #iteration.best memory
      individual.ib <- which.max(sapply(bf.results, function(x) return(x$solution.phe$pheromone)))
      logged.ib <- bf.results[[individual.ib]]
      solution.ib <- list()
      for (i in seq_along(short.factor.structure)) {
        solution.ib[[i]] <- seq_along(short.factor.structure[[i]])%in%bf.results[[individual.ib]]$selected[[i]]
      }
      phe.ib <- bf.results[[individual.ib]]$solution.phe$pheromone
      selected.ib <- bf.results[[individual.ib]]$selected
      
      #updated global best
      if (inherits(objective, 'stuartEmpiricalObjective')) {
        if (run > max(c(burnin, 1))) {
          if(all(is.na(logged.gb$solution.phe[-1]))) phe.gb <- 0
          else phe.gb <- do.call(objective$func, logged.gb$solution.phe[-1])
        }
      }
      
      #global.best memory
      if (phe.ib > phe.gb | generation == 1) {
        solution.gb <- solution.ib
        logged.gb <- logged.ib
        phe.gb <- phe.ib
        selected.gb <- selected.ib
      }
      
      #create matrix of combinations already evaluated
      tried <- rbind(tried, combi_mat)
      
      #create log
      log <- c(log, lapply(bf.results, function(x) c(run = run, x)))
      utils::setTxtProgressBar(progress,generation)

      # check for convergence
      
      if (generation >= generations_cur) {
        end.reason <- 'Maximum number of generations exceeded.'
        break
      }
      
      reinit <- FALSE
      conv <- FALSE
      qual.ib <- c(qual.ib, phe.ib)
      
      if ('variance' %in% reinit.criterion) {
        if (generation > max(min(c(.1*generations_cur, 10)),1) & stats::var(qual.ib/qual.ib[1]) <= reinit.tolerance_va_cur) {
          reinit <- TRUE
        }
      }
      
      if ('variance' %in% convergence.criterion) {
        if (generation > max(min(c(.1*generations_cur, 10)),1) & stats::var(qual.ib/qual.ib[1]) <= tolerance_va_cur) {
          conv <- TRUE
        }
      }
      
      if ('median' %in% reinit.criterion) {
        if (generation > max(min(c(.1*generations_cur, 10)),1) & ((phe.ib - stats::median(pheromones))/phe.ib) <= reinit.tolerance_md_cur) {
          reinit <- TRUE
        }
      }
      
      if ('median' %in% convergence.criterion) {
        if (generation > max(min(c(.1*generations_cur, 10)),1) & ((phe.ib - stats::median(pheromones))/phe.ib) <= tolerance_md_cur) {
          conv <- TRUE
        }
      }
      
      if ('geno.within' %in% c(convergence.criterion, reinit.criterion)) {
        geno_var <- list()
        tmp <- c(0, cumsum(unlist(capacity)))
        for (i in 1:length(short.factor.structure)) {
          all_items <- seq_along(short.factor.structure[[i]])
          sel_items <- combi_mat[, ((tmp[i]+1):tmp[i+1])]
          geno[[i]] <- t(apply(sel_items, 1, function(x) all_items %in% x))
          geno_var[[i]] <- colMeans(geno[[i]])*(1-colMeans(geno[[i]]))
        }
        if ('geno.within' %in% convergence.criterion &
            all(unlist(sapply(geno_var, function(x) x < (tolerance_gw_cur*(1-tolerance_gw_cur)))))) {
          conv <- TRUE
        }
        if ('geno.within' %in% reinit.criterion &
            all(unlist(sapply(geno_var, function(x) x < (reinit.tolerance_gw_cur*(1-reinit.tolerance_gw_cur)))))) {
          reinit <- TRUE
        }
      }
      
      if ('geno.between' %in% c(convergence.criterion, reinit.criterion)) {
        if (run == 1) {
          geno_r1 <- vector('list', length(short.factor.structure))
          geno_r1 <- lapply(geno_r1, function(x) 0)
          geno_d1 <- geno_r1
        } else {
          geno_r1 <- geno
        }
        geno_d2 <- geno_d1
        tmp <- c(0, cumsum(unlist(capacity)))
        for (i in 1:length(short.factor.structure)) {
          all_items <- seq_along(short.factor.structure[[i]])
          sel_items <- combi_mat[, ((tmp[i]+1):tmp[i+1])]
          geno[[i]] <- colMeans(t(apply(sel_items, 1, function(x) all_items %in% x)))
          geno_d1[[i]] <- abs(geno_r1[[i]]-geno[[i]])
        }
        
        if ('geno.between' %in% convergence.criterion &
            all(sapply(geno_d1, function(x) all(x < tolerance_gb_cur))) &
            all(sapply(geno_d2, function(x) all(x < tolerance_gb_cur)))) {
          conv <- TRUE
        }
        if ('geno.between' %in% reinit.criterion &
            all(sapply(geno_d1, function(x) all(x < reinit.tolerance_gb_cur))) &
            all(sapply(geno_d2, function(x) all(x < reinit.tolerance_gb_cur)))) {
          reinit <- TRUE
        }
      }
      
      # update empirical objective
      if (inherits(objective, 'stuartEmpiricalObjective') & run > burnin) {
        args <- c(objective$call, x = list(log))
        objective <- do.call(empiricalobjective, args)
      }
      
      
      # reinitialization
      if (reinit & cur_reinit.n > 0) {
        keep <- max(round(individuals_cur * (1 - reinit.prop_cur)), round(individuals_cur * elitism_cur))
        n <- individuals_cur - keep
        combinations <- do.call('generate.combinations',mget(names(formals(generate.combinations))))
        for (i in 1:length(short.factor.structure)) {
          nextgen[[i]] <- rbind(nextgen[[i]][1:keep, ], combinations$combi[[i]])
        }
        cur_reinit.n <- cur_reinit.n - 1
        message('\nReinitialized population. Generation counter reset.')
        generation <- 1
        utils::setTxtProgressBar(progress,generation)
      } else { # convergence
        if (conv) {
          end.reason <- 'Algorithm converged.'
          break
        } else {
          generation <- generation + 1
        }
      }
      
      #check for duplication in nextgen
      duplicate <- match(data.frame(t(do.call(cbind, nextgen))), data.frame(t(tried)))
      filter <- data.frame(matrix(which(is.na(duplicate)),
        nrow=sum(is.na(duplicate)),
        ncol=length(short.factor.structure)))
      
      #go on to next generation
      combi <- nextgen
      run <- run + 1

    } # end loop
     
    #feedback
    close(progress)
    message(paste('\nSearch ended.',end.reason))      
    
    # reformat log
    #generate matrix output
    mat_fil <- c('lvcor', 'lambda', 'theta', 'psi', 'alpha', 'beta', 'nu')
    mat_fil <- mat_fil[mat_fil %in% names(formals(objective$func))]
    mats <- as.list(vector('numeric', length(mat_fil)))
    names(mats) <- mat_fil
    
    for (m in seq_along(mat_fil)) {
      mats[[m]] <- sapply(log, function(x) x$solution.phe[mat_fil[m]])
      names(mats[[m]]) <- 1:length(log)
    }
    
    # apply final pheromone function retroactively (empirical objectives)
    if (inherits(objective, 'stuartEmpiricalObjective')) {
      final_pheromone <- sapply(log, function(x) {
        if (x$solution.phe$pheromone == 0) 0
        else {do.call(objective$func, x$solution.phe[-1])}
      })
    }

    tmp <- sapply(log, `[[`, 1)
    tmp <- unlist(lapply(table(tmp), function(x) seq(1, x)))
    log <- cbind(cumsum(tmp == 1), tmp, t(sapply(log, function(x) array(data=unlist(x$solution.phe[!names(x$solution.phe)%in%mat_fil])))))
    log <- data.frame(log)
    names(log) <- c('run', 'ind',names(bf.results[[1]]$solution.phe)[!names(bf.results[[1]]$solution.phe)%in%mat_fil])
    if (inherits(objective, 'stuartEmpiricalObjective')) {
      log$pheromone <- final_pheromone
    }
    
    #return to previous random seeds
    if (!is.null(seed)) {
      RNGkind(old.kind)
      .Random.seed <<- old.seed
    }
    
    # Preparing output
    convergence <- vector('list', length(convergence.criterion))
    names(convergence) <- convergence.criterion
    if ('variance' %in% names(convergence)) convergence[['variance']] <- stats::var(qual.ib/qual.ib[1])
    if ('median' %in% names(convergence)) convergence[['median']] <- phe.ib - stats::median(pheromones)
    if ('geno.within' %in% names(convergence)) convergence[['geno.within']] <- lapply(geno, colMeans)
    if ('geno.between' %in% names(convergence) & 'geno_d1'%in%ls()) convergence[['geno.between']] <- geno_d1
    
    tolerance <- list(variance = tolerance_va, median = tolerance_md, 
      geno.within = tolerance_gw, geno.between = tolerance_gb)
    reinit.tolerance <- list(variance = reinit.tolerance_va, median = reinit.tolerance_md, 
      geno.within = reinit.tolerance_gw, geno.between = reinit.tolerance_gb)
    
    #Generating Output
    tmp <- c(0, cumsum(unlist(capacity)))
    for (i in 1:length(short.factor.structure)) {
      all_items <- seq_along(short.factor.structure[[i]])
      sel_items <- combi_mat[, ((tmp[i]+1):tmp[i+1])]
      geno[[i]] <- colMeans(t(apply(sel_items, 1, function(x) all_items %in% x)))
    }
    names(geno) <- names(short.factor.structure)
    
    for (i in seq_along(solution.gb)) names(solution.gb[[i]]) <- short.factor.structure[[i]]
    names(solution.gb) <- names(short.factor.structure)
    
    results <- mget(grep('.gb',ls(),value=TRUE))
    results$selected.items <- translate.selection(selected.gb,factor.structure,short)
    results$log <- log
    results$log_mat <- mats
    results$pheromones <- pheromones
    results$parameters <- list(generations, individuals, selection, selection.pressure, elitism, mutation, mating.index, mating.size,
      mating.criterion, immigration, convergence.criterion, tolerance, 
      reinit.n, reinit.criterion,
      reinit.tolerance, reinit.prop,
      seed, objective, factor.structure)
    names(results$parameters) <- c('generations', 'individuals', 'selection', 'selection.pressure', 'elitism', 'mutation', 
      'mating.index', 'mating.size', 'mating.criterion', 'immigration', 'convergence.criterion', 'tolerance', 
      'reinit.n', 'reinit.criterion', 'reinit.tolerance', 'reinit.prop',
      'seed', 'objective', 'factor.structure')
    results$convergence <- convergence
    results$genotype <- geno
    results$end.reason <- end.reason

    return(results)
    
  }