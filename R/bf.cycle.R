bf.cycle <-
function(run,
  filter, combi,
  data, auxi,                                                   #data and selection coding
  capacity,
  long.equal, comparisons.equal, comparisons.invariance,
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
    run.options <- mget(names(formals(paste('run',software,sep='.'))))
    
    if (length(comparisons.equal)==0) {
      run.options$output.model <- FALSE
      solution.fit <- do.call(paste('run',software,sep='.'), run.options)
    } else {
      run.options$output.model <- TRUE
      solution.fit <- do.call(paste('run',software,sep='.'), run.options)

      if (!is.na(solution.fit)[[1]]) {
        comps <- list()
        for (i in seq_along(comparisons.equal)) {
          run.options$long.equal <- comparisons.equal[[i]]
          run.options[grep('invariance',names(run.options))] <- comparisons.invariance[[i]]
          comp.fit <- do.call(paste('run',software,sep='.'), run.options)
          comps <- c(comps, unlist(compute.comparisons(objective, comp.fit, solution.fit, names(comparisons.equal)[i])))
          if (is.logical(all.equal(objective.preset.comparisons, objective))) {
            names(comps)[grepl('delta\\.', names(comps))] <- gsub(paste0('\\.', names(comparisons.equal)[[i]]), '', names(comps)[grepl('delta\\.', names(comps))])
          }
        }
        solution.fit <- solution.fit[names(solution.fit)!='model']
        solution.fit <- c(solution.fit, unlist(comps))
      }
    }
  }

  #compute pheromone
  fitness.options <- as.list(formals(fitness))
  fitness.options <- mget(names(fitness.options))
  if (any(sapply(mtmm,length))>1) fitness.options$criteria <- c(as.character(fitness.options$criteria)[-1],'con')
  solution.phe <- do.call(fitness,fitness.options)
  if (!is.null(objective)) {
    if ('rel'%in%names(formals(objective))) {
      if (all(is.na(solution.phe$rel))) solution.phe$rel <- rep(NA,length(factor.structure)*max(c(1,sum(!is.na(unique(data[, grouping]))))))
    }
  }
  
  return(list(selected=selected,solution.phe=solution.phe))

}
