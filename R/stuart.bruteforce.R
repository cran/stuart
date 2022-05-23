stuart.bruteforce <-
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

  objective=NULL, ignore.errors=FALSE,                        #fitness function
  
  suppress.model=FALSE, analysis.options=NULL,                   #Additional modeling
  
  filename,

  ...                                                            #All the other stuff
) { #start function

  log <- NULL

  #Give Feedback about combinations
  message('Generating all possible combinations.')

  # Generate Combinations
  full <- TRUE
  n <- NULL
  combinations <- do.call('generate.combinations',mget(names(formals(generate.combinations))))

  filter <- combinations$filter
  combi <- combinations$combi

  #creating user feedback
  message('Running STUART with Brute-Force.\n')
  progress <- utils::txtProgressBar(style=3)
  utils::setTxtProgressBar(progress,0)
  count.gb <- 0

  output.model <- FALSE
  svalues <- FALSE
  bf.args <- mget(names(formals(bf.cycle))[-1])
  
  if (software=='Mplus') {
    #file location
    if (is.null(filename)) filename <- paste0(tempdir(), '/stuart')
    
    #writing the data file
    utils::write.table(data,paste(filename,'_data.dat',sep=''),
      col.names=FALSE,row.names=FALSE,na='-9999',
      sep='\t',dec='.')
  }
  
  #parallel processing for R-internal estimations
  if (software=='lavaan') {
    if (cores>1) {
      if (.Platform$GUI=='RStudio') message('\nProgressbars are not functional when utilizing multiple cores for bruteforce in RStudio.')
      #set up parallel processing on windows
      if (grepl('Windows',Sys.info()[1],ignore.case=TRUE)) {
        if (!.Platform$GUI=='RStudio') message('\nProgressbars are not functional when utilizing multiple cores for bruteforce in Windows.')
        cl <- parallel::makeCluster(cores)
        
        bf.results <- parallel::parLapply(cl,1:nrow(filter),function(run) {
          utils::setTxtProgressBar(progress, ceiling(run/(10*cores))/(nrow(filter)/(10*cores)));
          do.call('bf.cycle',c(run,bf.args))
        })
        parallel::stopCluster(cl)
      }
      
      #run ants in parallel on unixies
      else {
        bf.results <- parallel::mclapply(1:nrow(filter),
          function(run) {     
            utils::setTxtProgressBar(progress, ceiling(run/(10*cores))/(nrow(filter)/(10*cores)));
              do.call('bf.cycle',c(run,bf.args))
          },
          mc.cores=cores
        )
      }
    }
    
    else {
      bf.results <- lapply(1:nrow(filter),
        function(run) {     
          utils::setTxtProgressBar(progress, run/nrow(filter));
          do.call('bf.cycle',c(run,bf.args))
        }
      )
    }
  }
  
  #serial processing if Mplus is used (Mplus-internal parallelization is used)
  if (software=='Mplus') {
    bf.args$filename <- filename
    bf.args$cores <- cores
    bf.results <- lapply(1:nrow(filter),
      function(run) {     
        utils::setTxtProgressBar(progress, run/nrow(filter));
        do.call('bf.cycle',c(run,bf.args))
      }
    )
  }
  
  # Evaluate using empirical objective
  if (inherits(objective, 'stuartEmpiricalObjective')) {
    args <- c(objective$call, x = list(bf.results))
    objective <- do.call(empiricalobjective, args)
    bf.results <- lapply(bf.results, function(x) {
      x$solution.phe$pheromone <- do.call(objective$func, x$solution.phe[-1])
      return(x)})
  }
  
    #generate matrix output
  mat_fil <- c('lvcor', 'lambda', 'theta', 'psi', 'alpha', 'beta', 'nu')
  mat_fil <- mat_fil[mat_fil %in% names(formals(objective$func))]
  mats <- as.list(vector('numeric', length(mat_fil)))
  names(mats) <- mat_fil
  
  for (m in seq_along(mat_fil)) {
    mats[[m]] <- sapply(bf.results, function(x) x$solution.phe[mat_fil[m]])
    names(mats[[m]]) <- 1:nrow(filter)
  }
  
  log <- cbind(1:nrow(filter),t(sapply(bf.results, function(x) array(data=unlist(x$solution.phe[!names(x$solution.phe)%in%mat_fil])))))
  log <- data.frame(log)
  names(log) <- c('run',names(bf.results[[1]]$solution.phe)[!names(bf.results[[1]]$solution.phe)%in%mat_fil])

  #best solution
  run.gb <- which.max(sapply(bf.results, function(x) return(x$solution.phe$pheromone)))
  if (length(run.gb) > 1) {
    warning('The highest pheromone was achieved by multiple solutions. Only the first is reported.',call.=FALSE)
    run.gb <- run.gb[1]
  }
  phe.gb <- bf.results[[run.gb]]$solution.phe$pheromone
  selected.gb <- bf.results[[run.gb]]$selected

  close(progress)
  message('\nSearch ended.')

  tried <- try(do.call(cbind, lapply(filter, function(y) do.call(rbind,lapply(y, function(x) combi[[1]][x, ])))), silent = TRUE)
  if (inherits(tried, 'try-error')) warning('The full list of evaluated solutions could not be retrieved.',call.=FALSE)

  # construction solution in standard format
  solution.gb <- short.factor.structure
  for (i in 1:length(short.factor.structure)) {
    solution.gb[[i]] <- seq_along(short.factor.structure[[i]])
    solution.gb[[i]] <- solution.gb[[i]] %in% selected.gb[[i]]
    names(solution.gb[[i]]) <- short.factor.structure[[i]]
  }
  
  results <- mget(grep('.gb',ls(),value=TRUE))
  results$selected.items <- translate.selection(selected.gb,factor.structure,short)
  results$log <- log
  results$log_mat <- mats
  results$tried <- tried
  results$pheromones <- NULL
  results$parameters <- list(objective=objective, factor.structure=factor.structure)
  return(results)

}
