data.prep <-
function(
  data, factor.structure,                                               #simple prerequisites
  
  capacity=NULL, item.weights=NULL, item.invariance='congeneric',       #items
  repeated.measures=NULL, long.invariance='strict',                     #longitudinal relations
  mtmm=NULL, mtmm.invariance='configural',                              #mtmm relations
  grouping=NULL, group.invariance='strict',                             #grouping relations
  comparisons=NULL,
  auxiliary=NULL,                                                       #add variables
  objective=NULL,
  ...
) { #function begin

  #check for correct variable naming
  tmp <- c(unlist(factor.structure),grouping) #all names provided

  if (any(!is.element(tmp,names(data)))) {
    stop(paste0('The variable names you provided do not match the variable names in your dataset. Missing variables are: ', paste(tmp[which(!tmp%in%names(data))], sep = ', ') ,'.\n'),call.=FALSE)
  }

  #create phantom longitudinal data, if only cross-sectional
  if (is.null(repeated.measures)) {
    repeated.measures <- as.list(names(factor.structure))
    names(repeated.measures) <- names(factor.structure)
    long.invariance <- 'configural'
  } else {
    tmp <- as.list(c(repeated.measures,setdiff(names(factor.structure),unlist(repeated.measures))))
    names(tmp) <- c(names(repeated.measures),setdiff(names(factor.structure),unlist(repeated.measures)))
    repeated.measures <- tmp
    if (length(long.invariance)!=1 & length(long.invariance)!=length(repeated.measures)) {
      stop('The number of longitudinal invariance levels and the number of factors are not compatible.\n',call.=FALSE)
    }
  }

  #create phantom mtmm data, if only one method
  if (is.null(mtmm)) {
    mtmm <- as.list(names(factor.structure))
    names(mtmm) <- names(factor.structure)
    mtmm.invariance <- 'configural'
  } else {
    tmp <- as.list(c(mtmm,setdiff(names(factor.structure),unlist(mtmm))))
    names(tmp) <- c(names(mtmm),setdiff(names(factor.structure),unlist(mtmm)))
    mtmm <- tmp
    if (length(mtmm.invariance)!=1 & length(mtmm.invariance)!=length(mtmm)) {
      stop('The number of MTMM invariance levels and the number of factors are not compatible.\n',call.=FALSE)
    }
  }

  #create a longitudinal factor structure
  long.factor.structure <- factor.structure
  long.factor.structure[!(names(long.factor.structure)%in%sapply(repeated.measures,function(x) x[1]))] <- NULL

  #create an mtmm factor structure
  mtmm.factor.structure <- factor.structure
  if (mtmm.invariance!='none') {
    mtmm.factor.structure[!(names(mtmm.factor.structure)%in%sapply(mtmm,function(x) x[1]))] <- NULL
    
    #create a short factor structure (minimal)
    short.factor.structure <- factor.structure[intersect(names(long.factor.structure),names(mtmm.factor.structure))]
  } else {
    short.factor.structure <- long.factor.structure
  }
  

  #create a list of short allocation
  short <- as.list(names(short.factor.structure))
  names(short) <- unlist(short)
  
  for (i in 1:length(short.factor.structure)) {
    counter <- 1
    repeat {
      filter <- sapply(repeated.measures,function(x)  x[any(short[[i]]%in%x)])
      short[[i]] <-union(short[[i]],unlist(filter))
      if (mtmm.invariance!='none') {
        filter <- sapply(mtmm,function(x)  x[any(short[[i]]%in%x)])
        short[[i]] <-union(short[[i]],unlist(filter))
      }
      if (counter-length(short[[i]])==0) break
      counter <- length(short[[i]])
    }
  }
  
  #### MOVE? ####
  #create a vector of items
  items <- unlist(factor.structure,use.names=FALSE)
  
  #create the subset of auxiliary variables
  auxi <- data[,auxiliary, drop = FALSE]
  names(auxi) <- auxiliary
  ####       ####

  #create vectors of invariance assumptions
  item.invariance <- as.list(array(item.invariance,length(short.factor.structure)))
  item.invariance <- lapply(item.invariance, function(x) ordered(x, levels = c('congeneric', 'ess.equivalent', 'equivalent', 'ess.parallel', 'parallel')))

  inv.ordering <- function(x) ordered(x, levels = c('none', 'configural', 'weak', 'strong', 'strict'))
  long.invariance <- as.list(array(long.invariance,length(long.factor.structure)))
  long.invariance <- lapply(long.invariance, inv.ordering)
  mtmm.invariance <- as.list(array(mtmm.invariance,length(mtmm.factor.structure)))
  mtmm.invariance <- lapply(mtmm.invariance, inv.ordering)
  group.invariance <- lapply(group.invariance, inv.ordering)
  
  #errors for wrong invariance settings
  if (any(is.na(c(unlist(long.invariance), unlist(mtmm.invariance), unlist(group.invariance))))) {
    stop(paste0('Invariance levels across repeated measurements, groups, and sources of information must be one of ',paste(levels(long.invariance[[1]]),collapse=', '),'.'),call.=FALSE)
  }
  
  if (any(is.na(unlist(item.invariance)))) {
    stop(paste0('Item invariance must be one of ',paste(levels(item.invariance[[1]]),collapse=', '),'.'),call.=FALSE)
  }
  
  # Expanding the capacity
  if (is.numeric(capacity)) {
    capacity <- as.list(rep(capacity,length(short.factor.structure)))
  } else {
    capacity <- capacity[names(factor.structure)%in%names(short.factor.structure)]
  }
  
  # #create a vector of group item invariance assumptions
  # if (length(item.group.invariance)!=1 & length(item.group.invariance)!=length(factor.structure)) {
  #   stop('The number of group invariance levels and the number of factors are not compatible.\n',call.=FALSE)
  # }
  # 
  # item.group.invariance <- as.list(array(item.group.invariance,length(mtmm.factor.structure)))
  number.of.some <- lapply(short.factor.structure,length)
  label.change <- FALSE
  inv.args <- mget(names(formals(invariance.implementation)))
  
  #implement invariances of items
  long.equal <- do.call(invariance.implementation, inv.args)
  #implement invariances for comparison models
  comparisons.equal <- list()
  comparisons.invariance <- list()
  for (i in comparisons) {
    cur.inv <- mget(paste0(i,'.invariance'))
    tmp.args <- inv.args
    if (i == 'item') {
      if (all(unlist(tmp.args[[paste0(i,'.invariance')]])=='congeneric')) {
        stop('The assumed item invariance is congeneric for all facets. Invariance testing via comparisons is not possible in this case.', call.=FALSE)
      }
      tmp.args[[paste0(i,'.invariance')]] <- lapply(tmp.args[[paste0(i,'.invariance')]], function(x) ordered(levels(x)[max(as.numeric(x)-1, 1)], levels = c('congeneric', 'ess.equivalent', 'equivalent', 'ess.parallel', 'parallel')))
    } else {
      tmp.args[[paste0(i,'.invariance')]] <- lapply(tmp.args[[paste0(i,'.invariance')]], function(x) {
        if (as.numeric(x)>2) inv.ordering(levels(x)[as.numeric(x)-1])
        else x
      })
    }
    comparisons.equal[[i]] <- do.call(invariance.implementation, tmp.args)
    comparisons.invariance[[i]] <- tmp.args[grep('invariance', names(tmp.args))]
    
  }
  
  #set preset for objective functions
  if (is.null(objective)) {
    if (is.null(comparisons)) objective <- fixedobjective()
    else objective <- fixedobjective(comparisons = comparisons)
  }
  if (inherits(objective, 'function')) {
    objective <- list(func = objective, string = toString(body(objective)[-1]))
    class(objective) <- 'stuartManualObjective'
  }
  
  output <- list(short.factor.structure,short,long.equal,comparisons.equal,comparisons.invariance,
      capacity,data,factor.structure,auxi,item.invariance,
      repeated.measures,long.invariance,mtmm,mtmm.invariance,grouping,group.invariance,
      objective)
  names(output) <- c('short.factor.structure','short','long.equal','comparisons.equal','comparisons.invariance',
      'capacity','data','factor.structure','auxi','item.invariance',
      'repeated.measures','long.invariance','mtmm','mtmm.invariance','grouping','group.invariance',
      'objective')

  return(output)

}
