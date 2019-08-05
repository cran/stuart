data.prep <-
function(
  data, factor.structure,                                               #simple prerequisites
  
  capacity=NULL, item.weights=NULL, item.invariance='congeneric',       #items
  repeated.measures=NULL, long.invariance='strict',                     #longitudinal relations
  mtmm=NULL, mtmm.invariance='configural',                              #mtmm relations
  grouping=NULL, group.invariance='strict',                             #grouping relations

  auxiliary=NULL,                                                       #add variables
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
  long.invariance <- as.list(array(long.invariance,length(long.factor.structure)))
  mtmm.invariance <- as.list(array(mtmm.invariance,length(mtmm.factor.structure)))
  
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
  
  #implement invariances of items
  long.equal <- invariance.implementation(data,
    factor.structure,short.factor.structure,short,
    long.factor.structure,repeated.measures,
    mtmm.factor.structure,mtmm,
    lapply(short.factor.structure,length),
    item.invariance,long.invariance,mtmm.invariance,
    group.invariance,
    grouping)
  
  output <- list(short.factor.structure,short,long.equal,
      capacity,data,factor.structure,auxi,item.invariance,
      repeated.measures,long.invariance,mtmm,mtmm.invariance,grouping,group.invariance)
  names(output) <- c('short.factor.structure','short','long.equal',
      'capacity','data','factor.structure','auxi','item.invariance',
      'repeated.measures','long.invariance','mtmm','mtmm.invariance','grouping','group.invariance')

  return(output)

}
