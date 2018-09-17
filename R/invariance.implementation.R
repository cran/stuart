invariance.implementation <-
function(
  data,                                                         #data
  factor.structure, short.factor.structure, short,
  long.factor.structure, repeated.measures,                     #data structure
  mtmm.factor.structure, mtmm,
  
  number.of.some,

  invariance, long.invariance, mtmm.invariance, group.invariance, #invariances

  grouping,                                                     #additional data
  label.change=FALSE                                            #replace label names?
) { #begin function

  #errors for wrong invariance settings
  tmp.inv <- c('none','configural','weak','strong','strict')
  if (any(c((!unlist(long.invariance)%in%tmp.inv),(!unlist(mtmm.invariance)%in%tmp.inv),(!unlist(group.invariance)%in%tmp.inv)))) {
    stop(paste0('Invariance levels across repeated measurements, groups, and sources of information must be one of ',paste(tmp.inv,collapse=', '),'.'),call.=FALSE)
  }
  
  tmp.inv <- c('congeneric','ess.equivalent','equivalent','ess.parallel','parallel')
  if (any(!unlist(invariance%in%tmp.inv))) {
    stop(paste0('Item invariance must be one of ',paste(tmp.inv,collapse=', '),'.'),call.=FALSE)
  }
  
  #invariance parameters
  equal <- rep(list(list(lam=NA,alp=NA,eps=NA)),length(factor.structure))
  names(equal) <- names(factor.structure)
  
  
  for (i in 1:length(factor.structure)) {
    locate <- which(unlist(lapply(short,
      function(x) is.element(names(factor.structure)[i],x))))
    
    locate.long <- which(unlist(lapply(repeated.measures,
      function(x) is.element(names(factor.structure)[i],x))))
    locate.mtmm <- which(unlist(lapply(mtmm,
      function(x) is.element(names(factor.structure)[i],x))))    

    locate.long2 <- which(repeated.measures[[which(unlist(lapply(repeated.measures,
      function(x) is.element(names(factor.structure)[i],x))))]]==names(factor.structure)[[i]])
    locate.mtmm2 <- which(mtmm[[which(unlist(lapply(mtmm,
      function(x) is.element(names(factor.structure)[i],x))))]]==names(factor.structure)[[i]])
    
    
    if (is.element(names(factor.structure)[i],names(short.factor.structure))) {
      
      #generate indices (item, construct, method, occasion)
      equal[[i]]  <- lapply(equal[[i]],function(x) array(NA,c(number.of.some[[locate]],4)))
      
      #item/subtest indices
      if (invariance[[locate]]=='congeneric') {
        equal[[i]]$lam[,1] <- 1:nrow(equal[[i]]$lam)
      } else {
        equal[[i]]$lam[,1] <- 1
      }
      
      if (invariance[[locate]]%in%c('equivalent','parallel')) {
        equal[[i]]$alp[,1] <- 1
      } else {            
        equal[[i]]$alp[,1] <- 1:nrow(equal[[i]]$alp)
      }
      
      if (invariance[[locate]]%in%c('ess.parallel','parallel')) {
        equal[[i]]$eps[,1] <- 1
      } else {
        equal[[i]]$eps[,1] <- 1:nrow(equal[[i]]$eps)
      }
    
      #construct indices
      equal[[i]]$lam[,2] <- locate
      equal[[i]]$alp[,2] <- locate
      equal[[i]]$eps[,2] <- locate
      
      #method indices
      equal[[i]]$lam[,3] <- locate.mtmm2
      equal[[i]]$alp[,3] <- locate.mtmm2
      equal[[i]]$eps[,3] <- locate.mtmm2
      
      #occasion indices
      equal[[i]]$lam[,4] <- locate.long2
      equal[[i]]$alp[,4] <- locate.long2
      equal[[i]]$eps[,4] <- locate.long2
      
    } else {
      equal[[i]] <- equal[[names(locate)]]
      
      if (mtmm.invariance[[locate.mtmm]]!='strict') equal[[i]]$eps[,3] <- locate.mtmm2
      
      if (mtmm.invariance[[locate.mtmm]]%in%c('weak','configural')) equal[[i]]$alp[,3] <- locate.mtmm2
      
      if (mtmm.invariance[[locate.mtmm]]=='configural') equal[[i]]$lam[,3] <- locate.mtmm2

      if (long.invariance[[locate.long]]!='strict') equal[[i]]$eps[,4] <- locate.long2

      if (long.invariance[[locate.long]]%in%c('weak','configural')) equal[[i]]$alp[,4] <- locate.long2
      
      if (long.invariance[[locate.long]]=='configural') equal[[i]]$lam[,4] <- locate.long2
      
    }
  }

  for (i in 1:length(equal)) {
    for (j in 1:length(equal[[i]])) {
      equal[[i]][[j]] <- paste0(names(equal[[i]])[j],apply(equal[[i]][[j]],1,paste0,collapse=''))
    }
  }
    

  #implementing group invariance
  if (!is.null(grouping)) {
    group <- as.factor(data[,grouping])
    group <- droplevels(group)
    equal <- list(equal,equal)

    for (i in 2:length(levels(group))) {
      equal[[i]] <- equal[[1]] }

    #add variable residuals
    if (group.invariance!='strict') {
      for (i in 2:length(levels(group))) {
        for (j in 1:length(factor.structure)) {
          equal[[i]][[j]]$eps <- paste(equal[[i]][[j]]$eps,'g',i,sep='')
        }
      }
    }
    
    #add variable intercepts
    if (group.invariance%in%c('weak','configural')) {
      for (i in 2:length(levels(group))) {
        for (j in 1:length(factor.structure)) {
          equal[[i]][[j]]$alp <- paste(equal[[i]][[j]]$alp,'g',i,sep='')
        }
      }
    }

    #add variable loadings
    if (group.invariance=='configural') {
      for (i in 2:length(levels(group))) {
        for (j in 1:length(factor.structure)) {
          equal[[i]][[j]]$lam <- paste(equal[[i]][[j]]$lam,'g',i,sep='')
        }
      }
    }
  }
  
  if (label.change) {
    tmp <- utils::as.relistable(equal)
    tmp <- unlist(tmp)
    tmp <- gsub('lam','gam',tmp)
    tmp <- gsub('alp','mu',tmp)
    tmp <- gsub('eps','zet',tmp)
    equal <- utils::relist(tmp)
  }

  return(equal)
}