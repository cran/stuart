bf.cycle <-
function(run,
  filter, combi,
  data, auxi,                                                   #data and selection coding
  capacity,
  long.equal,
  factor.structure, repeated.measures, mtmm, grouping,          #basic requirements
  short.factor.structure, short,
  item.invariance, long.invariance, mtmm.invariance, group.invariance, #invariance settings
  analysis.options, suppress.model,                             #additional analysis options
  objective,                                                 #fitness function to call
  software,output.model=FALSE,svalues=FALSE,
  ignore.errors=FALSE,
  filename,cores,
  ...
) {#function begin

  selected <- list()
  for (j in 1:length(combi)) {
    selected[[j]] <- as.numeric(combi[[j]][filter[run,j],])
  }
  selected.items <- translate.selection(selected,factor.structure,short)
    
  if (any(duplicated(unlist(selected.items)))) {
    solution.fit <- NA
  } else {
    output.model=FALSE
    run.options <- names(formals(paste('run',software,sep='.')))
    solution.fit <- do.call(paste('run',software,sep='.'),as.list(mget(run.options)))
  }

  #compute pheromone
  fitness.options <- as.list(formals(fitness))
  fitness.options <- mget(names(fitness.options))
  if (any(sapply(mtmm,length))>1) fitness.options$criteria <- c(as.character(fitness.options$criteria)[-1],'con')
  solution.phe <- do.call(fitness,fitness.options)
  if (!is.null(objective)) {
    if ('rel'%in%names(formals(objective))) {
      if (all(is.na(solution.phe$rel))) solution.phe$rel <- rep(NA,length(factor.structure)*max(c(1,length(unique(data[,grouping])))))
    }
  }
  
  return(list(selected=selected,solution.phe=solution.phe))

}
