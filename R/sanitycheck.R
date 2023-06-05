sanitycheck <- function(data, factor.structure,capacity,
  repeated.measures,mtmm,grouping=NULL,
  objective=NULL,localization,software='lavaan',
  comparisons) {
  
  #sanity check
  if (any(sapply(data[, unlist(factor.structure)], function(x) all(class(x)=='factor')))) {
    if (!all(sapply(data[, unlist(factor.structure)], nlevels) %in% c(0, 2))) {
      stop('Currently only binary, ordinal, and continuous items are supported.', call. = FALSE)
    }
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
  
  if (!is.null(objective) & !inherits(objective, 'function')) {
    if (!is.function(objective$func)) {
      stop('The objective function you requested is not a function.',call.=FALSE)
    }
    if (is.null(mtmm)&'con'%in%names(formals(objective$func))) {
      stop('The objective function you requested uses consistency in the pheromone computation but mtmm=NULL',call.=FALSE)
    }
    if (any(sapply(data[, unlist(factor.structure)], is.factor)) & any(names(formals(objective$func))=='srmr') & software == 'Mplus') {
      stop('Mplus does not provide estimates for the SRMR when handling ordinal variables. Please change your objective function.', call. = FALSE)
    }
    
    if (any(sapply(data[, unlist(factor.structure)], is.factor)) & !any(grepl('.scaled|.robust', names(formals(objective$func)))) & software == 'lavaan') {
      warning('It is highly recommended to used either scaled or robust versions of model fit criteria in your objective function when modeling ordinal indicators with lavaan.', call. = FALSE)
    }
    
    if (any(grepl('delta.', names(formals(objective$func)))) & is.null(grouping) & is.null(mtmm) & is.null(repeated.measures)) {
      stop('Your objective function contains a model comparison value, but you provided no groups, mwthods, or occasions to compare.', call. = FALSE)
    }
  }
  
  if (!localization%in%c('arcs','nodes')) stop('Pheromones must be localized to arcs or nodes.',call.=FALSE)
  
  tmp_filt <- c(is.null(repeated.measures) & any(comparisons == 'long'), is.null(mtmm) & any(comparisons == 'mtmm'), is.null(grouping) & any(comparisons == 'group'))
  tmp_pot <- cbind(c('long', 'mtmm', 'group'), c('repeated.measures', 'mtmm', 'grouping'))
  if (any(tmp_filt)) {
    stop(paste0('You requested comparisons for \"', tmp_pot[tmp_filt, 1], '\" but did not provide all necessary information in \"', tmp_pot[tmp_filt, 2], '\".'), call. = FALSE)
  }
  
}