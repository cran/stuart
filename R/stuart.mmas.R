stuart.mmas <-
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
  
  ants=16, colonies=256, evaporation=.95,                        #general ACO parameters
  deposit='ib', pbest=.005, localization='nodes',                #MMAS parameters
  alpha=1, beta=1, pheromones=NULL, heuristics=NULL,
  tolerance=.5, schedule='run',                                  #tolerance for convergence

  suppress.model=FALSE, analysis.options=NULL,                   #Additional modeling
  seed,
  filename,

  ...                                                            #All the other stuff
) { #function begin

  #initialize fitness results of best solutions
  phe.ib <- 0
  phe.gb <- 0
  
  #initialize upper and lower limits
  phe.max <- 0
  phe.min <- 0
  
  #counting
  run <- 1
  colony <- 1

  #compute number of decisions and avg for limits
  deci <- sum(unlist(capacity))
  avg <- mean(c(sapply(short.factor.structure,length),(1+sapply(short.factor.structure,length)-unlist(capacity))))

  #recode deposit to numeric
  deposit_save <- deposit
  deposit[deposit=='ib'] <- 1
  deposit[deposit=='gb'] <- 2
  class(deposit) <- 'numeric'
  
  #initialize scheduling 
  scheduled <- c('ants','colonies','evaporation','pbest','alpha','beta','tolerance','deposit')
  
  #global assignment to avoid check note
  ants_cur <- NA
  colonies_cur <- NA
  evaporation_cur <- NA
  pbest_cur <- NA
  alpha_cur <- NA
  beta_cur <- NA
  tolerance_cur <- NA
  deposit_cur <- NA
  
  filt <- sapply(mget(scheduled),is.array)
  for (i in 1:length(scheduled[!filt])) {
    assign(paste0(scheduled[!filt][i],'_cur'),mget(scheduled[!filt][i])[[1]])
  }
  if (length(scheduled[filt])>0) {
    scheduled <- scheduled[filt]
    for (i in 1:length(scheduled)) {
      tmp <- mget(scheduled[i])[[1]]
      if (!any(c(0,1)%in%tmp[,1])) {
        stop(paste('The parameter schedule for',scheduled[i],'does not contain a value for the first colony.'),call.=FALSE)
      }
      tmp <- tmp[which.min(tmp[,1]),2]
      assign(paste0(scheduled[i],'_cur'),tmp)
      assign(scheduled[i],cbind(mget(scheduled[i])[[1]],FALSE))
    }
  } else {
    scheduled <- NULL
  }
  
  #initialize pheromones
  if (is.null(pheromones)) pheromones <- init.pheromones(short.factor.structure,localization,alpha_cur)

  if (is.null(heuristics)) {
    heuristics <- lapply(pheromones,function(x) x^(1/1e+100))
  } else {
    if (attr(heuristics,'localization')!=localization) {
      stop('The localization of the heuristics does not match your setting for localization.',call.=FALSE)
    }
  }

  
  #set random seed, if provided
  if (!is.null(seed)) {
    old.seed <- .Random.seed
    old.kind <- RNGkind()[1]
    RNGkind("L'Ecuyer-CMRG")
    set.seed(seed)
  }
  
  
  ### Loops ###
  log <- list()
  tried <- matrix(NA, ncol = sum(unlist(capacity)))[-1,]
  counter <- matrix(NA, ncol=2)[-1,]
  

  #creating user feedback
  message('Running STUART with MMAS.\n')
  progress <- utils::txtProgressBar(0,max(c(colonies,1)),style=3)
  count.gb <- 0

  repeat { #over colonies

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
        if (schedule=='colony') {
          if (any(tmp[,1]==colony)) {
            message(paste0('Scheduled value of ',scheduled[i],' updated to ',tmp[which(tmp[,1]==run),2],'.'))
          }
          tmp <- tmp[max(which(tmp[,1]<=colony)),2]
        }
        if (schedule=='mixed') {
          if (any(tmp[,3]-(tmp[,1]<=colony)<0)) mix_new <- TRUE
          if (mix_new) {
            tmp[,3] <- tmp[,3]|(tmp[,1]<=colony)
          } else {
            tmp[,3] <- tmp[,3]|(tmp[,1]<=(colony+max(c(tmp[as.logical(tmp[,3]),1],1)-1)))
          }
          assign(scheduled[i],tmp)
          if (any(tmp[,1]==colony)&mix_new) {
            message(paste0('Scheduled value of ',scheduled[i],' updated to ',tmp[which.max(tmp[as.logical(tmp[,3]),1]),2],'.'))
          }
          tmp <- tmp[which.max(tmp[as.logical(tmp[,3]),1]),2]
        }
        assign(paste0(scheduled[i],'_cur'),tmp)
      }
    }

    output.model <- FALSE
    svalues <- FALSE
    cons.args <- mget(names(formals(paste('construction',localization,sep='.'))))
    if (length(scheduled[scheduled%in%names(cons.args)])>0) {
      ant.args[scheduled[scheduled%in%names(cons.args)]] <- mget(paste(scheduled[scheduled%in%names(cons.args)],'cur',sep='_'))
    }
    constructed <- lapply(1:ants_cur, function(x) do.call(paste('construction',localization,sep='.'),cons.args))
    
    combi <- lapply(constructed, function(x) x$selected)
    combi <- do.call(Map, c(rbind, combi))
    duplicate <- match(data.frame(t(do.call(cbind, combi))), data.frame(t(tried)))
    filter <- data.frame(matrix(which(is.na(duplicate)),
      nrow=sum(is.na(duplicate)),
      ncol=length(short.factor.structure)))
    
    tried <- rbind(tried, do.call(cbind, combi))
    
    ant.args <- mget(names(formals(bf.cycle))[-1])
    if (length(scheduled[scheduled%in%names(ant.args)])>0) {
      ant.args[scheduled[scheduled%in%names(ant.args)]] <- mget(paste(scheduled[scheduled%in%names(ant.args)],'cur',sep='_'))
    }
    
    if (nrow(filter) > 0) {
      #parallel processing for R-internal estimations
      if (software=='lavaan') {
        if (cores>1) {
          #set up parallel processing on windows
          if (grepl('Windows',Sys.info()[1],ignore.case=TRUE)) {
            cl <- parallel::makeCluster(cores)
            if (!is.null(seed)) {
              seed_cur <- sample(1e+9,1)
              parallel::clusterSetRNGStream(cl,seed_cur)
            }
            ant.results <- parallel::parLapply(cl,1:nrow(filter),function(run) do.call(bf.cycle,c(run,ant.args)))
            parallel::stopCluster(cl)
          }
          
          #run ants in parallel on unixies
          else {
            ant.results <- parallel::mclapply(1:nrow(filter),function(run) do.call(bf.cycle,c(run,ant.args)), mc.cores=cores)
          }
        }
        
        #serial processing for single cores
        else {
          ant.results <- lapply(1:nrow(filter),function(run) do.call(bf.cycle,c(run,ant.args)))
        }
      }
      
      #serial processing if Mplus is used (Mplus-internal parallelization is used)
      if (software=='Mplus') {
        ant.args$filename <- filename
        ant.args$cores <- cores
        ant.results <- lapply(1:nrow(filter),function(run) do.call(bf.cycle,c(run,ant.args)))
      }
    }
    
    #fill in results for duplicates
    tmp <- vector('list', ants_cur)
    tmp[filter[,1]] <- ant.results
    tmp[sapply(tmp,is.null)] <- log[stats::na.omit(duplicate)]
    ant.results <- tmp
    
    #iteration.best memory
    ant.ib <- which.max(sapply(ant.results, function(x) return(x$solution.phe$pheromone)))
    solution.ib <- constructed[[ant.ib]]$solution
    phe.ib <- ant.results[[ant.ib]]$solution.phe$pheromone
    selected.ib <- ant.results[[ant.ib]]$selected

    #feedback
    utils::setTxtProgressBar(progress,colony)

    #global.best memory
    if (phe.ib > phe.gb | run == 1) {
      count.gb <- count.gb + 1
      solution.gb <- solution.ib
      phe.gb <- phe.ib
      selected.gb <- selected.ib

      # in cases of mixed counting, reset mix_new
      mix_new <- FALSE
      
      #new solution user feedback
      message(paste('\nGlobal best no.',count.gb,'found. Colony counter reset.'))
      
      #restart the count
      colony <- 1
      utils::setTxtProgressBar(progress,0)
    }

    else {
      colony <- colony + 1
    }

    #compute upper and lower limits
    phe.max <- phe.gb/(1-evaporation_cur)
    phe.min <- (phe.max*(1-pbest_cur^(1/deci)))/((avg-1)*pbest_cur^(1/deci))

    if (phe.min >= phe.max) {
      stop('The lower pheromone limit is larger than the upper pheromone limit. This may indicate that none of the initial solutions were viable due to estimation problems.\n',call.=FALSE)
    }
    
    #updated pheromones
    pheromones <- mmas.update(pheromones,phe.min,phe.max,evaporation_cur,localization,
      get(paste('phe',c('ib','gb')[deposit_cur],sep='.')),get(paste('solution',c('ib','gb')[deposit_cur],sep='.')))

    
    #create log
    log <- c(log, ant.results)
    counter <- rbind(counter, c(run, ants_cur))
    # log <- rbind(log,cbind(rep(run,ants_cur),1:ants_cur,t(sapply(ant.results, function(x) array(data=unlist(x$solution.phe))))))
    
    
    #check for convergence
    if (localization=='arcs') {
      conv <- lapply(pheromones,function(x) x[lower.tri(x)])
    } else {
      conv <- pheromones
    }
    tmp.min <- sapply(conv, function(x) sum(phe.min-tolerance_cur < x & x < phe.min+tolerance_cur))
    tmp.max <- sapply(conv, function(x) sum(phe.max-tolerance_cur < x & x < phe.max+tolerance_cur))
    tmp.all <- sapply(conv, length)
    tmp <- cbind(tmp.min,tmp.max)
    abort.sequence <- all(rowSums(tmp)==tmp.all) & all(tmp!=0)&run>1


    #abort if converged
    if (abort.sequence) {
      end.reason <- 'Algorithm converged.'
      break
    }
    if (colony > colonies_cur) {
      end.reason <- 'Maximum number of colonies exceeded.'
      break
    }

    #keep on counting
    run <- run + 1
  }

  #feedback
  close(progress)
  message(paste('\nSearch ended.',end.reason))      

  #return to previous random seeds
  if (!is.null(seed)) {
    RNGkind(old.kind)
    .Random.seed <<- old.seed
  }
  
  # reformat log
  tmp <- as.numeric(unlist(apply(counter, 1, function(x) seq(1,x[2]))))
  log <- cbind(cumsum(tmp==1),tmp,t(sapply(log, function(x) array(data=unlist(x$solution.phe)))))
  log <- data.frame(log)
  names(log) <- c('run','ant',names(ant.results[[1]]$solution.phe))
  
  
  #Generating Output
  results <- mget(grep('.gb',ls(),value=TRUE))
  results$selected.items <- translate.selection(selected.gb,factor.structure,short)
  results$log <- log
  results$pheromones <- pheromones
  results$parameters <- list(ants=ants,colonies=colonies,evaporation=evaporation,
    deposit=deposit_save,pbest=pbest,localization=localization,
    alpha=alpha,beta=beta,tolerance=tolerance,schedule=schedule,phe.max=phe.max,phe.min=phe.min,
    seed=seed,
    objective=objective,
    heuristics=heuristics,
    factor.structure=factor.structure)
  return(results)

}
