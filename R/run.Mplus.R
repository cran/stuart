run.Mplus <-
function(
  data, auxi, 
  capacity,
  selected, selected.items,
  long.equal,
  factor.structure, repeated.measures, grouping,
  short.factor.structure, short, mtmm,
  item.invariance, long.invariance, mtmm.invariance, group.invariance,
  
  analysis.options=NULL, suppress.model=FALSE,
  
  output.model=FALSE, svalues=FALSE,
  ignore.errors=FALSE,
  filename=NULL, cores
) { #begin function
  
  #prepare data for model fit
  model.data <- data[,unlist(selected.items)]
  if (!is.null(grouping)) {
    model.data$group <- data[,grouping]
    if (!is.numeric(data[, grouping])) {
      stop('Mplus requires numerical variables. Your grouping variable is not numeric.', call. = FALSE)
    }
  }
  model.data <- data.frame(model.data,auxi)

  #check for ordinal items
  ordinal <- any(sapply(data[, unlist(factor.structure)], function(x) class(x)[1]) == 'ordered')
  
  #file location
  if (is.null(filename)) filename <- paste0(tempdir(), '/stuart')
  
  # Check for inappropriate analysis options
  if (!is.null(analysis.options)) {
    tmp <- c('title', 'data', 'variable', 'analysis', 'model', 'constraints', 'output', 'savedata')
    names(analysis.options) <- tolower(names(analysis.options))
    if (any(!names(analysis.options) %in% tmp)) stop(paste('The analysis options supplied to Mplus must be a list named according to the Mplus input sections.'), call. = FALSE)
  }
  
  #define empty input
  input <- NULL

  #write Mplus "Title" section
  input <- paste0('Title: Subtest Construction using STUART \n', analysis.options$title, '\n')

  #write Mplus "Data" section
  input <- paste0(input,'Data: file=',filename,'_data.dat; \n', analysis.options$data, '\n')
  
  # Add Analysis options to "Data" Section
  input <- paste0(input, '\n', analysis.options$data, '\n')
  
  #write Mplus "Variable" section
  input <- paste0(input,'Variable: \n\tnames=')
  
  input <- paste0(input,paste(names(data),collapse=c('\n\t\t')),';\n')
  
  input <- paste0(input,'\tmissing=ALL(-9999);\n',
    '\tusevariables=')
  
  input <- paste0(input,paste(unlist(selected.items),collapse=c('\n\t\t')),';\n')
  
  if (any(sapply(data[, unlist(selected.items)], is.factor))) {
    input <- paste0(input,'\tcategorical=', 
      paste(unlist(selected.items)[sapply(data[, unlist(selected.items)], is.factor)],collapse=c('\n\t\t')),';\n')
  }
  
  if (!is.null(grouping)) {
    input <- paste(input,paste0('grouping = ', grouping, ' (',paste(stats::na.omit(unique(model.data$group)),stats::na.omit(unique(model.data$group)),sep='=',collapse=' '),')'),';\n')
  }
  
  input <- paste0(input,unlist(analysis.options[grepl('^vari*',names(analysis.options),ignore.case=TRUE)][1]),'\n')
  
  
  #write Mplus "Analysis" section
  input <- paste0(input,'Analysis: \n\tprocessors=',cores,';\n')
  
  input <- paste0(input, analysis.options$analysis,'\n')
  
  #write Mplus "Model" section
  input <- paste0(input,'Model:\n')
  
  #define mplus model syntax
  if (!suppress.model) {
    
    #write Mplus input (no groups)
    if (is.null(grouping)) {
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
        tmp.inv <- lapply(long.equal[[i]],function(x) return(x[tmp.sel]))
        if (ordinal) {
          tmp.inv$alp <- array(dim=0)
          for (l in tmp.sel) { #across selected items
            tmp.inv$alp <- c(tmp.inv$alp, long.equal[[i]]$alp[cthresh[l]:(cthresh[l+1]-1)])
          }
        }

        #factor loadings
        input <- paste(input,'\n',
                       names(selected.items[i]),'by',
                       paste0(tmp.sit,' (',tmp.inv$lam,')',collapse='\n\t\t'),';\n')
        
        #residual variances & intercepts
        tmp.thr <- c(0, cumsum(nthresh[tmp.sit]))
        for (j in seq_along(tmp.sit)) {
          if (is.factor(data[, tmp.sit[j]])) { #categorical
            input <- paste(input, 
              paste0('[', tmp.sit[j], '$', 1:nthresh[tmp.sit[j]], '] (', tmp.inv$alp[(tmp.thr[j] + 1):tmp.thr[j+1]], ');', sep = '', collapse = '\n'), sep = '\n')
          } else {
            input <- paste(input, 
              paste0('[', tmp.sit[j], '] (', tmp.inv$alp[(tmp.thr[j] + 1):tmp.thr[j+1]], ');', sep = '', collapse = '\n'), sep = '\n')
            input <- paste(input,
              paste0(tmp.sit[j],' (',tmp.inv$eps[j],');',collapse='\n'),sep='\n')
          }
        }
        
        #estimate latent regressions (MTMM)
        if (names(selected.items[i])%in%lapply(mtmm, function(x) x[1])) {
          tmp <- mtmm[[which(unlist(lapply(mtmm, function(x) x[1]))%in%names(selected.items[i]))]][-1]
          if (length(tmp)>0) {
            regs <- expand.grid(sapply(tmp,function(x) names(selected.items[x])),names(selected.items[i]))
            regs <- sapply(regs,as.character)
          
            if (is.null(nrow(regs))) {
              tmp <- paste0(paste(regs,collapse=' on '),';\n')
            } else {
              tmp <- paste0(paste(apply(regs,1,paste,collapse=' on '),collapse=';\n'),';\n')
            }
            
            input <- paste(input,tmp,sep='\n')
          }
        }
      }
      
      #set latent means
      for (i in 1:length(factor.structure)) {
        if (long.invariance[[which(unlist(lapply(repeated.measures,function(x) is.element(names(factor.structure)[i],x))))]]%in%c('strong','strict')) {
          if (names(selected.items[i])%in%lapply(repeated.measures, function(x) x[1])) {
            input <- paste(input,
              paste0('[',names(selected.items[i]),'@0];',collapse='\n'),sep='\n')
          } else {
            input <- paste(input,
              paste0('[',names(selected.items[i]),'*];',collapse='\n'),sep='\n')
          }
        }
      }
    }
    
    #write Mplus input (grouping)
    else {
      #write the (item) factor structure
      for (i in 1:length(selected.items)) { #over factors

        #shorten the writing by creating tmp-data
        tmp.fil <- which(unlist(lapply(short,
          function(x) is.element(names(factor.structure)[i],x))))
        tmp.sel <- selected[[tmp.fil]]
        tmp.sit <- selected.items[[i]]
        
        #factor loadings (overall)
        input <- paste(input,'\n',
                       names(selected.items[i]),'by',
                       paste0(tmp.sit,collapse='\n\t\t'),';\n')
        
        #residual variances & intercepts
        for (j in seq_along(tmp.sit)) {
          # for categorical
          if (is.factor(data[, tmp.sit[j]])) {
            input <- paste(input, 
              paste0('[',tmp.sit[j],'$1];', collapse = '\n'),sep='\n')
          } 
          # for continuous
          else {
            input <- paste(input,
              paste0(tmp.sit[j],';',collapse='\n'),sep='\n')
            input <- paste(input,
              paste0('[',tmp.sit[j],'];', collapse = '\n'),sep='\n')
          }
        }
        
        #estimate latent regressions (MTMM)
        if (names(selected.items[i])%in%lapply(mtmm, function(x) x[1])) {
          tmp <- mtmm[[which(unlist(lapply(mtmm, function(x) x[1]))%in%names(selected.items[i]))]][-1]
          if (length(tmp)>0) {
            regs <- expand.grid(sapply(tmp,function(x) names(selected.items[x])),names(selected.items[i]))
            regs <- sapply(regs,as.character)
            
            if (is.null(nrow(regs))) {
              tmp <- paste0(paste(regs,collapse=' on '),';\n')
            } else {
              tmp <- paste0(paste(apply(regs,1,paste,collapse=' on '),collapse=';\n'),';\n')
            }
            
            input <- paste(input,tmp,sep='\n')
          }
        }
      }

      #group specific models
      for (k in 1:length(long.equal)) { #over groups
        
        #write grouping header
        input <- paste(input,'\n',
                       'Model',stats::na.omit(unique(model.data$group))[k],':\n')
        
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
          
          tmp.thr <- c(0, cumsum(nthresh[tmp.sit]))

          tmp.inv <- list(NA)
          tmp.inv <- lapply(long.equal[[k]][[i]],function(x) return(x[tmp.sel]))
          if (ordinal) {
            tmp.inv$alp <- character()
            for (j in seq_along(tmp.sel)) {
              tmp.inv$alp <- c(tmp.inv$alp, long.equal[[k]][[i]]$alp[cthresh[tmp.sel[j]]:(cthresh[tmp.sel[j]]+(nthresh[tmp.sel[j]]-1))])
            }
          }

          #factor loadings
          tmp.lam <- paste(names(selected.items[i]),'by',
                           paste0(tmp.sit[1],'@1\n\t\t'))
          tmp.lam <- paste(tmp.lam,
                           paste0(tmp.sit[-1],' (',tmp.inv$lam[-1],')',collapse='\n\t\t'),';\n')
          input <- paste(input,'\n',tmp.lam)
          
          #residual variances & intercepts
          for (j in seq_along(tmp.sit)) {
            if (is.factor(data[, tmp.sit[j]])) { #categorical
              input <- paste(input, 
                paste0('[', tmp.sit[j], '$', 1:nthresh[tmp.sit[j]], '] (', tmp.inv$alp[(tmp.thr[j] + 1):tmp.thr[j+1]], ');', sep = '', collapse = '\n'), sep = '\n')
            } else {
              input <- paste(input, 
                paste0('[', tmp.sit[j], '] (', tmp.inv$alp[(tmp.thr[j] + 1):tmp.thr[j+1]], ');', sep = '', collapse = '\n'), sep = '\n')
              input <- paste(input,
                paste0(tmp.sit[j],' (',tmp.inv$eps[j],');',collapse='\n'),sep='\n')
            }
          }
        }

        #set latent means
        for (i in 1:length(factor.structure)) {
          if (long.invariance[[which(unlist(lapply(repeated.measures,function(x) is.element(names(factor.structure)[i],x))))]]%in%c('strong','strict')) {
            if (names(selected.items[i])%in%lapply(repeated.measures, function(x) x[1])&
                k==1) {
              input <- paste(input,
                paste0('[',names(selected.items[i]),'@0];',collapse='\n'),sep='\n')
            } else {
              input <- paste(input,
                paste0('[',names(selected.items[i]),'*];',collapse='\n'),sep='\n')
            }
          }
        }
      }    
    }
  }
  
  input <- paste0(input, analysis.options$model, '\n')
  
  # Add model constraint section of analysis options
  input <- paste0(input, analysis.options$constraints, '\n')
  
  #write Mplus "Output" section
  if (!output.model | ordinal) {
    input <- paste(input,'Output: STDYX Tech4 NOSERROR;\n')

    # Add Analysis Options Output section
    input <- paste0(input, analysis.options$output, '\n')
  }
  if (ordinal & output.model) warning('Due to a bug in MplusAutomation standard errors can currently not be computed when using ordinal indicators.', call. = FALSE)
  
  else {
    input <- paste(input,'Output: STDYX Tech4;\n')
    # Add Analysis Options Output section
    input <- paste0(input, analysis.options$output, '\n')
  }
  
  # write "Savedata" section
  input <- paste0(input, analysis.options$savedata, '\n')
    
  #create Mplus input file
  cat(input,file=paste0(filename,'.inp'))
  
  # Begin silencer (for noisy MplusAutomation)
  tmp <- textConnection(NULL, 'w')
  sink(tmp)
  
  # Run Model
  MplusAutomation::runModels(paste0(filename, '.inp'))
  
  #import Mplus output
  MplusOut <- try(suppressWarnings(MplusAutomation::readModels(paste0(filename, '.out'))), silent = TRUE)
  
  # end silencer
  sink()
  close(tmp) 

  if (class(MplusOut)[1] == 'try-error') return(output = list(NA))
  
  if (svalues) {
    tmp <- readLines(paste0(filename, '.out'))
    tmp <- tmp[(grep('USED AS STARTING VALUES',tmp)+1):(grep('^TECHNICAL',tmp)[1]-1)]
    MplusOut$svalues <- tmp
  }

  exclusion <- FALSE
  if (length(MplusOut$errors) > 0) return(output=list(NA))
  if (!ignore.errors) {
    exclusion <- any(sapply(MplusOut$warnings, function(x) any(grepl('POSITIVE DEFINITE|NO CONVERGENCE|CHECK YOUR MODEL', x))))
  }

  #return list of NA if errors occurred
  if (exclusion) {
    return(output=list(NA))
  } else {

    #extract the fit statistics reported by Mplus
    output <- as.list(MplusOut$summaries)
    
    name <- c('rmsea','srmr','wrmr','cfi','tli','chisq','df','pvalue','aic','bic','abic','npar')
    locator <- c('RMSEA_Estimate', 'SRMR', 'WRMR', 'CFI', 'TLI', 'ChiSqM_Value', 'ChiSqM_DF', 'ChiSqM_PValue', 'AIC', 'BIC', 'aBIC', 'Parameters')
    
    for (i in seq_along(locator)) {
      names(output)[which(names(output)==locator[i])] <- name[i]
    }
    
    #extract latent correlations
    if (is.null(grouping)) {
      lvcor <- list(MplusOut$tech4$latCorEst)
      psi <- list(MplusOut$tech4$latCovEst)
    } else {
      lvcor <- lapply(MplusOut$tech4, function(x) x$latCorEst)
      psi <- lapply(MplusOut$tech4, function(x) x$latCovEst)
    }

    # workaround for WLSMV bug in MplusAutomation
    if (!is.null(lvcor)) {
      lvcor <- lapply(lvcor, function(x) {
        x[upper.tri(x)] <- t(x)[upper.tri(x)]
        dimnames(x) <- list(names(selected.items), names(selected.items))
        return(x)})
      psi <- lapply(psi, function(x) {
        x[upper.tri(x)] <- t(x)[upper.tri(x)]
        dimnames(x) <- list(names(selected.items), names(selected.items))
        return(x)})
      names(psi) <- names(lvcor) <- NULL
    }

    output$lvcor <- lvcor
    
    # compute rho estimate of reliability
    esti <- MplusOut$parameters$unstandardized
    if (is.null(grouping)) esti$Group <- 1
    
    nitems <- length(unlist(selected.items))
    nfacto <- length(selected.items)
    
    # Extract lambda
    lambda <- lapply(unique(esti$Group), function(x) {
      tmp <- esti[grepl('.BY$', esti$paramHeader) & esti$Group == x, ]
      tmp$lv <- gsub('.BY$', '', tmp[grep('.BY$', tmp$paramHeader), 'paramHeader'])
      tmp <- stats::reshape(tmp[, c('param', 'est', 'lv')], timevar = 'lv', idvar = 'param', direction = 'wide')
      tmp <- tmp[, -1]
      tmp[is.na(tmp)] <- 0
      tmp <- as.matrix(tmp)
      dimnames(tmp) <- list(unlist(selected.items), names(selected.items))
      return(tmp)
    })
    
    # Extract alpha
    tmp <- esti[(grepl('Intercepts', esti$paramHeader) & esti$param %in% toupper(unlist(selected.items))) | grepl('Thresholds', esti$paramHeader), ]
    alpha <- tapply(tmp$est, tmp$Group, identity, simplify = FALSE)

    # Extract theta
    var_types <- sapply(model.data, function(x) class(x)[1])
    var_types <- var_types[names(var_types)!='group']
    
    if (all(var_types%in%c('numeric','integer'))) {
      tmp <- esti[grepl('Residual.Variances', esti$paramHeader) & esti$param %in% toupper(unlist(selected.items)), ]
      theta <- tapply(tmp$est, tmp$Group, function(x) diag(x), simplify = FALSE)
    } else {
      tmp <- MplusOut$parameters$r2
      if (is.null(grouping)) {
        tmp$Group <- 1
      }
      theta <- tapply(tmp$resid_var, tmp$Group, function(x) diag(x), simplify = FALSE)
      
      if (any(var_types%in%c('numeric','integer'))) {
        for (i in seq_along(theta)) {
          tmp <- esti[grepl('Residual.Variances', esti$paramHeader) & esti$param %in% toupper(unlist(selected.items)), ]
          tmp <- tmp[tmp$Group == unique(tmp$Group)[i], 'est']
          theta[[i]][is.na(theta[[i]])] <- tmp
        } 
      }
    }
    
    dimnames(alpha) <- dimnames(theta) <- names(lambda) <- NULL

    # Compute Reliabilities and Composite Reliabilities    
    rel <- lapply(lambda,function(x) rep(NA,ncol(x)))
    crel <- rep(NA,length(lambda))
    for (i in 1:length(rel)) {
      for (j in 1:length(rel[[i]])) {
        filter <- which(lambda[[i]][,j]!=0)
        rel[[i]][j] <- sum(lambda[[i]][,j,drop=FALSE]%*%psi[[i]][j,j,drop=FALSE]%*%t(lambda[[i]][,j,drop=FALSE]))/(sum(lambda[[i]][,j,drop=FALSE]%*%psi[[i]][j,j,drop=FALSE]%*%t(lambda[[i]][,j,drop=FALSE]))+sum(theta[[i]][filter,filter,drop=FALSE]))
      }
      # workaround for absence of short.factor.structure when crossvalidating
      if (class(try(short.factor.structure,silent=TRUE))=='try-error') {
        short.factor.structure <- as.list(rep(NA,ncol(lambda[[i]])))
        names(short.factor.structure) <- colnames(lambda[[i]])
      }
      reffilter <- colnames(lambda[[i]])%in%names(short.factor.structure)
      filter <- rowSums(lambda[[i]][,reffilter,drop=FALSE]!=0)>0
      
      crel[i] <- sum(lambda[[i]][,reffilter,drop=FALSE]%*%psi[[i]][reffilter,reffilter,drop=FALSE]%*%t(lambda[[i]][,reffilter,drop=FALSE]))/(sum(lambda[[i]][,reffilter,drop=FALSE]%*%psi[[i]][reffilter,reffilter,drop=FALSE]%*%t(lambda[[i]][,reffilter,drop=FALSE]))+sum(theta[[i]][filter,filter,drop=FALSE]))
    }
    crel <- mean(crel)
    output$crel <- crel
    if (is.null(grouping)) {
      output$lambda <- lambda[[1]]
      output$theta <- theta[[1]]
      output$psi <- psi[[1]]
      output$alpha <- alpha[[1]]
      output$rel <- rel[[1]]
    } else {
      output$lambda <- lambda
      output$theta <- theta
      output$psi <- psi
      output$alpha <- alpha
      output$rel <- rel
    }
    
    tmp <- MplusOut$parameters$r2
    if (is.null(grouping)) tmp$Group <- 1
    tmp <- tmp[tmp$param %in% toupper(names(selected.items)), ]
    con <- mean(tapply(tmp$est, tmp$Group, mean))

    if (output.model) output$model <- MplusOut
    
    return(output=output)
  }  
  
} #end function