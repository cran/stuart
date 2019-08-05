sanitycheck <- function(data, factor.structure,capacity,
  repeated.measures,mtmm,
  objective=NULL,localization,software='lavaan') {
  
  if (is.null(objective)) objective <- objective.preset
  #sanity check
  if (any(sapply(data[, unlist(factor.structure)], function(x) all(class(x)=='factor')))) {
    if (!all(sapply(data[, unlist(factor.structure)], nlevels) %in% c(0, 2))) {
      stop('Currently only binary, ordinal, and continuous items are supported.', call. = FALSE)
    }
  }
  
  if (any(sapply(data[, unlist(factor.structure)], is.factor)) & any(names(formals(objective))=='srmr') & software == 'Mplus') {
    stop('Mplus does not provide estimates for the SRMR when handling ordinal variables. Please change your objective function.', call. = FALSE)
  }
  
  if (any(duplicated(names(factor.structure)))) {
    stop('You have provided duplicates in the name of factor.structure.',call.=FALSE)
  }
  
  if (is.null(capacity) | (is.numeric(capacity)&length(capacity)!=1) | (is.list(capacity)&length(capacity)!=length(factor.structure))) {
    stop('You must provide either one global value or a list of the same length as the factor.structure for capacity.',call.=FALSE)
  } 

  if (any(!unlist(repeated.measures)%in%names(factor.structure))) {
    stop(paste('One or more factors appearing in repeated.measures is not present in the factor.structure:',
      paste(unlist(repeated.measures)[!unlist(repeated.measures)%in%names(factor.structure)],collapse=', ')),call.=FALSE)
  }
  
  if (any(!unlist(mtmm)%in%names(factor.structure))) {
    stop(paste('One or more factors appearing in mtmm is not present in the factor.structure:',
      unlist(mtmm)[!unlist(mtmm)%in%names(factor.structure)]),call.=FALSE)
  }
  
  if (!is.null(objective)) {
    if (!is.function(objective)) {
      stop('The objective function you requested is not a function.',call.=FALSE)
    }
    if (is.null(mtmm)&'con'%in%names(formals(objective))) {
      stop('The objective function you requested uses consistency in the pheromone computation but mtmm=NULL',call.=FALSE)
    }
  }
  
  if (!localization%in%c('arcs','nodes')) stop('Pheromones must be localized to arcs or nodes.',call.=FALSE)
  
}