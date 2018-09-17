ant.cycle <-
function(
  localization='arcs',                                          #type of solution construction
  data, auxi, use.order, pheromones,                            #data and selection coding
  alpha, beta, heuristics,
  capacity,
  long.equal,                                                   #invariance labels
  factor.structure, repeated.measures, mtmm, grouping,          #basic requirements
  short.factor.structure, short,
  item.invariance, long.invariance, mtmm.invariance, group.invariance, #invariance settings
  analysis.options, suppress.model,                             #additional analysis options
  objective,                                                 #fitness function to call
  software,output.model=FALSE,svalues=FALSE,
  ignore.errors=FALSE,
  filename,cores
) { #begin function

  .Deprecated('bf.cycle')
  
  #construct solution
  tmp <- mget(names(formals(paste('construction',localization,sep='.'))))
  constructed <- do.call(paste('construction',localization,sep='.'),tmp)

  solution <- constructed$solution
  selected <- constructed$selected
  
  #translate selection to names
  tmp <- mget(names(formals(translate.selection)))
  selected.items <- do.call(translate.selection,tmp)

  #specify modeling options
  run.options <- names(formals(paste('run',software,sep='.')))
  solution.fit <- do.call(paste('run',software,sep='.'),as.list(mget(run.options)))
  
  #compute pheromone
  fitness.options <- as.list(formals(fitness))
  fitness.options$solution.fit <- solution.fit
  fitness.options$objective <- objective
  if (any(sapply(mtmm,length))>1) fitness.options$criteria <- c(as.character(fitness.options$criteria)[-1],'con')
  solution.phe <- do.call(fitness,fitness.options)
  if (!is.null(objective)) {
    if ('rel'%in%names(formals(objective))) {
      if (all(is.na(solution.phe$rel))) solution.phe$rel <- rep(NA,length(factor.structure)*max(c(1,length(unique(data[,grouping])))))
    }
  }

  return(list(solution=solution,selected=selected,solution.phe=solution.phe))

}
