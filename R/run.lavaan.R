run.lavaan <-
function(
  data, auxi, 
  capacity,
  selected, selected.items,
  long.equal,
  factor.structure, repeated.measures, grouping,
  short.factor.structure, short, mtmm = NULL,
  item.invariance, long.invariance, mtmm.invariance, group.invariance,

  analysis.options = NULL, suppress.model = FALSE,

  output.model = FALSE,
  objective = NULL, ignore.errors = FALSE
) { #begin function

  #prepare data for model fit
  model.data <- data[,unlist(selected.items)]
  model.data$group <- data[,grouping]
  model.data <- data.frame(model.data,auxi)

  #define empty lavaan input
  input <- NULL

  #check for ordinal items
  ordinal <- any(sapply(data[, unlist(factor.structure)], function(x) class(x)[1]) == 'ordered')
  
  if (!suppress.model) {
    
    #write the (item) factor structure
    for (i in 1:length(selected.items)) { #over factors
      #shorten the writing by creating tmp-data
      tmp.fil <- which(unlist(lapply(short,
        function(x) is.element(names(factor.structure)[i],x))))
      tmp.sel <- selected[[tmp.fil]]
      tmp.sit <- selected.items[[i]]

      locate <- which(unlist(lapply(short,
        function(x) is.element(names(factor.structure)[i],x))))
      nthresh <- sapply(data[, factor.structure[[i]]], function(x) max(nlevels(x), 2))-1
      cthresh <- c(0, cumsum(nthresh))+1

      #write the labels (no grouping)
      if (is.null(grouping)) {
        tmp.inv <- lapply(long.equal[[i]],function(x) return(x[tmp.sel]))
        if (ordinal) {
          tmp.inv$alp <- array(dim=0)
          for (l in tmp.sel) { #across selected items
            tmp.inv$alp <- c(tmp.inv$alp, long.equal[[i]]$alp[cthresh[l]:(cthresh[l+1]-1)])
          }
        }
      }

      #write the labels (grouping)
      else {
        tmp.inv <- list(NA)
        for (k in 1:length(long.equal)) { #over groups
          tmp.inv[[k]] <- lapply(long.equal[[k]][[i]],function(x) return(x[tmp.sel]))
          if (ordinal) {
            tmp.inv[[k]]$alp <- character()
            for (j in seq_along(tmp.sel)) {
              tmp.inv[[k]]$alp <- c(tmp.inv[[k]]$alp, long.equal[[k]][[i]]$alp[cthresh[tmp.sel[j]]:(cthresh[tmp.sel[j]]+(nthresh[tmp.sel[j]]-1))])
            }
          }
          tmp.inv[[k]] <- unlist(tmp.inv[[k]])
        }
        tmp.inv <- data.frame(lapply(tmp.inv,data.frame))
        tmp.inv <- apply(tmp.inv,1,paste,collapse=',')
        tmp.inv <- paste('c(',tmp.inv,')',sep='')
        tmp.inv <- list(lam=grep('lam', tmp.inv, value = TRUE),
          alp=grep('alp', tmp.inv, value = TRUE),
          eps=grep('eps', tmp.inv, value = TRUE))
      }

      #factor loadings
      input <- paste(input,'\n',
        names(selected.items[i]),'=~',
        paste(tmp.inv$lam,'*',tmp.sit,sep='',collapse=' + '),sep='')

      #residual variances
      tmp_filt <- sapply(data[, tmp.sit], is.numeric)
      if (any(tmp_filt)) {
        input <- paste0(input, '\n',
          paste0(tmp.sit[tmp_filt], '~~', tmp.inv$eps[tmp_filt], '*', tmp.sit[tmp_filt], sep = '', collapse = '\n'))
      }

      #intercepts
      tmp.thr <- c(0, cumsum(nthresh[tmp.sit]))
      for (j in seq_along(tmp.sit)) {
        if (is.factor(data[, tmp.sit[j]])) {
          input <- paste(input, 
            paste(tmp.sit[j], '|', tmp.inv$alp[(tmp.thr[j] + 1):tmp.thr[j+1]], '*t', 1:nthresh[tmp.sit[j]], sep = '', collapse = '\n'), sep = '\n')
        } else {
          input <- paste(input, 
            paste(tmp.sit[j], '~', tmp.inv$alp[(tmp.thr[j] + 1):tmp.thr[j+1]], '*1', sep = '', collapse = '\n'), sep = '\n')
        }
      }

        
      #supress correlations between traits and methods (for CTC(M-1) structure)
#       if (names(selected.items[i])%in%lapply(mtmm, function(x) x[1])) {
#         tmp <- mtmm[[which(unlist(lapply(mtmm, function(x) x[1]))%in%names(selected.items[i]))]][-1]
#         tmp <- outer(names(selected.items[[i]]),sapply(tmp,function(x) names(selected.items[[x]])),
#           paste,sep=' ~~ 0*')
#         tmp <- paste(tmp,collapse='\n')
#         input <- paste(input,tmp,sep='\n')
#       }
  
      #estimate latent regressions (MTMM)
      if (names(selected.items[i])%in%lapply(mtmm, function(x) x[1])) {
        tmp <- mtmm[[which(unlist(lapply(mtmm, function(x) x[1]))%in%names(selected.items[i]))]][-1]
        if (length(tmp) > 0) {
          regs <- expand.grid(sapply(tmp,function(x) names(selected.items[x])),names(selected.items[i]))
          regs <- sapply(regs,as.character)
          if (is.null(nrow(regs))) {
            tmp <- paste(regs,collapse='~')
          } else {
            tmp <- paste(apply(regs,1,paste,collapse='~'),collapse='\n')
          }
          
          input <- paste(input,tmp,sep='\n')
        }
      }
      
      # write latent means
      if (long.invariance[[which(unlist(lapply(repeated.measures,function(x) is.element(names(factor.structure)[i],x))))]]%in%c('strong','strict')) {
        if (names(selected.items[i])%in%lapply(repeated.measures, function(x) x[1])) {
          if (!is.null(grouping)&unlist(group.invariance)%in%c('strong','strict')) {
            input <- paste(input,
              paste0(names(selected.items[i]),'~c(', paste(c(0,rep(NA,nlevels(as.factor(model.data$group))-1)),collapse=','),')*1',collapse='\n'),sep='\n')
          } else {
            input <- paste(input,
              paste0(names(selected.items[i]),'~ 0*1',collapse='\n'),sep='\n')
          }
        } else {
          input <- paste(input,
            paste0(names(selected.items[i]),'~ 1;',collapse='\n'),sep='\n')
        }
      }
      
      if (!is.null(grouping) & length(repeated.measures[[locate]])==1) {
        if (unlist(group.invariance) %in% c('strong', 'strict')) {
          input <- paste(input,
            paste0(names(selected.items[i]),'~c(', paste(c(0,rep(NA,nlevels(as.factor(model.data$group))-1)),collapse=','),')*1',collapse='\n'),sep='\n')
        }
      }
    }
  }

  if (is.data.frame(analysis.options$model)) {
    if (suppress.model) input <- analysis.options$model
    else input <- rbind(lavaan::lavParTable(input),analysis.options$model)
  }
  else input <- paste(input,analysis.options$model,sep='\n')
  
  #list of arguments to pass to lavaan
  if (is.null(analysis.options)) {
    analysis.options <- list(NULL) 
  }
    
  analysis.options$model <- input
  analysis.options$data <- model.data

  if (!is.null(grouping)) {
    analysis.options$group <- 'group'
  }
  if (is.null(analysis.options$missing) & !any(sapply(data[, unlist(factor.structure)], is.factor))) {
    analysis.options$missing <- 'fiml'
  }

  if (!output.model) {
    analysis.options$se <- 'none'
  }
  
  if (all(names(formals(objective$func))%in%c('crel', 'rel', 'con', 'lvcor', 'alpha', 'beta', 'lambda', 'psi', 'theta'))) {
    analysis.options$h1 <- FALSE
  }
  
  #imply sem() presets
  presets <- list(int.ov.free=TRUE,int.lv.free=FALSE,auto.fix.first=TRUE,std.lv=FALSE,
    auto.fix.single=TRUE,auto.var=TRUE,auto.cov.lv.x=TRUE,auto.th=TRUE,auto.delta=TRUE,auto.cov.y=TRUE,parameterization='delta')
  for (i in names(presets)) {
    if (!i %in% names(analysis.options)) analysis.options[i] <- presets[i]
  }
  
  #retain only the options that are accepted by lavaan
  analysis.options <- analysis.options[!sapply(analysis.options,is.null)]

  # Temporary, for multiple lavaan Versions. Will be removed.  
  if (utils::packageVersion('lavaan')<'0.5.23') {
    analysis.options <- analysis.options[is.element(names(analysis.options),names(formals(lavaan::lavaan)))]
  } else {
    analysis.options <- analysis.options[is.element(names(analysis.options),c(names(formals(lavaan::lavaan)),names(lavaan::lavOptions())))]
  }
  
  # tmp.cfa <- get('cfa',asNamespace('lavaan'))  
  output <- try(suppressWarnings(do.call(lavaan::lavaan,analysis.options)),silent=TRUE)

  if (inherits(output, 'try-error')) {
    warning('The lavaan input generated an error.',call.=FALSE)
    return(output=list(NA))
  }

  if (inherits(output, 'lavaan')) {

    if (!ignore.errors) {
      if (!suppressWarnings(lavaan::inspect(output, 'post.check'))) return(output = list(NA))
    }

    fits <- names(formals(objective$func))
    fits <- gsub('delta\\.', '', fits)
    fits <- gsub('\\.items|\\.long|\\.mtmm|\\.group', '', fits)
    fits <- unique(fits)
    fit <- try(suppressWarnings(lavaan::fitMeasures(output, fits)),silent=TRUE)
    
    if (inherits(fit, 'try-error')) {
      return(output=list(NA))
      warning('The lavaan estimation generated an error, most likely non-convergence.')
    } else {
      # compute composite reliability (overall)
      if (is.null(grouping)) {
        tmp <- lavaan::inspect(output,'est')
        
        theta <- tmp$theta[rownames(tmp$theta)%in%unlist(selected.items),colnames(tmp$theta)%in%unlist(selected.items), drop = FALSE]
        psi <- tmp$psi[rownames(tmp$psi)%in%names(factor.structure),colnames(tmp$psi)%in%names(factor.structure), drop  = FALSE]
        lambda <- tmp$lambda[rownames(tmp$lambda)%in%unlist(selected.items),colnames(tmp$lambda)%in%names(factor.structure), drop = FALSE]
        
        rel <- rep(NA,ncol(lambda))        
        for (i in 1:ncol(lambda)) {
          filter <- which(lambda[,i]!=0)
          rel[i] <- sum(lambda[,i,drop=FALSE]%*%psi[i,i,drop=FALSE]%*%t(lambda[,i,drop=FALSE]))/(sum(lambda[,i,drop=FALSE]%*%psi[i,i,drop=FALSE]%*%t(lambda[,i,drop=FALSE]))+sum(theta[filter,filter,drop=FALSE]))
        }
        # workaround for absence of short.factor.structure when crossvalidating
        if (inherits(try(short.factor.structure,silent=TRUE), 'try-error')) {
          short.factor.structure <- as.list(rep(NA,ncol(lambda)))
          names(short.factor.structure) <- colnames(lambda)
        }
        reffilter <- colnames(lambda)%in%names(short.factor.structure)
        filter <- rowSums(lambda[,reffilter,drop=FALSE]!=0)>0
        
        crel <- sum(lambda[,reffilter,drop=FALSE]%*%psi[reffilter,reffilter,drop=FALSE]%*%t(lambda[,reffilter,drop=FALSE]))/(sum(lambda[,reffilter,drop=FALSE]%*%psi[reffilter,reffilter,drop=FALSE]%*%t(lambda[,reffilter,drop=FALSE]))+sum(theta[filter,filter,drop=FALSE]))
        
        # pass matrices from lavaan to output
        theta <- tmp$theta
        psi <- tmp$psi
        lambda <- tmp$lambda
        alpha <- tmp$alpha
        nu <- tmp$nu
        beta <- tmp$beta
        if (ordinal) {
          tau <- tmp$tau
          delta <- tmp$delta
        }
        
        tmp <- lavaan::inspect(output,'rsquare')
        con <- mean(tmp[names(tmp)%in%names(factor.structure)])
          
      } else {
        tmp <- lavaan::inspect(output,'est')
        theta <- lapply(tmp,function(x) x$theta[rownames(x$theta)%in%unlist(selected.items),colnames(x$theta)%in%unlist(selected.items), drop = FALSE])
        psi <- lapply(tmp,function(x) x$psi[rownames(x$psi)%in%names(factor.structure),colnames(x$psi)%in%names(factor.structure), drop  = FALSE])
        lambda <- lapply(tmp,function(x) x$lambda[rownames(x$lambda)%in%unlist(selected.items),colnames(x$lambda)%in%names(factor.structure), drop = FALSE])
        
        rel <- lapply(lambda,function(x) rep(NA,ncol(x)))
        crel <- rep(NA,length(lambda))
        for (i in 1:length(rel)) {
          for (j in 1:length(rel[[i]])) {
            filter <- which(lambda[[i]][,j]!=0)
            rel[[i]][j] <- sum(lambda[[i]][,j,drop=FALSE]%*%psi[[i]][j,j,drop=FALSE]%*%t(lambda[[i]][,j,drop=FALSE]))/(sum(lambda[[i]][,j,drop=FALSE]%*%psi[[i]][j,j,drop=FALSE]%*%t(lambda[[i]][,j,drop=FALSE]))+sum(theta[[i]][filter,filter,drop=FALSE]))
          }
          # workaround for absence of short.factor.structure when crossvalidating
          if (inherits(try(short.factor.structure,silent=TRUE), 'try-error')) {
            short.factor.structure <- as.list(rep(NA,ncol(lambda[[i]])))
            names(short.factor.structure) <- colnames(lambda[[i]])
          }
          reffilter <- colnames(lambda[[i]])%in%names(short.factor.structure)
          filter <- rowSums(lambda[[i]][,reffilter,drop=FALSE]!=0)>0
          
          crel[i] <- sum(lambda[[i]][,reffilter,drop=FALSE]%*%psi[[i]][reffilter,reffilter,drop=FALSE]%*%t(lambda[[i]][,reffilter,drop=FALSE]))/(sum(lambda[[i]][,reffilter,drop=FALSE]%*%psi[[i]][reffilter,reffilter,drop=FALSE]%*%t(lambda[[i]][,reffilter,drop=FALSE]))+sum(theta[[i]][filter,filter,drop=FALSE]))
        }
        crel <- mean(crel)

        # pass matrices from lavaan to output
        alpha <- lapply(tmp,function(x) x$alpha)
        nu <- lapply(tmp,function(x) x$nu)
        beta <- lapply(tmp,function(x) x$beta)
        theta <- lapply(tmp,function(x) x$theta)
        lambda <- lapply(tmp,function(x) x$lambda)
        psi <- lapply(tmp,function(x) x$psi)
        if (ordinal) {
          tau <- lapply(tmp,function(x) x$tau)
          delta <- lapply(tmp,function(x) x$delta)
        }
        
        tmp <- lavaan::inspect(output,'rsquare')
        con <- mean(sapply(tmp,function(x) mean(x[!names(x)%in%names(model.data)])))
      }
      
      
      # Export the latent variable correlation matrix
      lvcor <- lavaan::inspect(output,'cor.lv')
      
      # Export optimizer information
      technicals <- lavaan::inspect(output, 'optim')
      technicals <- technicals[names(technicals) != 'control']

      #Listed in detail for quick overview of exported output
      results <- as.list(fit)
      results$crel <- crel
      results$rel <- rel
      results$lvcor <- lvcor
      results$lambda <- lambda
      results$theta <- theta
      results$psi <- psi
      results$alpha <- alpha
      results$nu <- nu
      results$beta <- beta
      if (!is.na(con)) results$con <- con
      if (ordinal) {
        results$tau <- tau
        results$delta <- delta
      }
      results$technicals <- technicals
      
      if (output.model) {
        results$model <- output
      }

      return(output=results)
    }
  }

}
