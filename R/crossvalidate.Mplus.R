### Function definition ----
crossvalidate.Mplus <-
function(
  selection, old.data, new.data, output.model=TRUE, 
  analysis.options = NULL, filename, ...
) { # begin function
  
  if (is.null(filename)) filename <- paste0(tempdir(), '/stuart')
  model <- selection$final$input$model
  model[grep('by\\s.+\\s\\(lam', model)] <- gsub(' \\(lam', '@1 \\(lam', model[grep('by\\s.+\\s\\(lam', model)])

  if (!is.null(selection$call$grouping)) {
    stop('Crossvalidate is currently not available for situations with multiple groups.', call. = FALSE)
  } else {
    grouping  <- 'stuart_sample'
    old.data[, grouping] <- 1
    new.data[, grouping] <- 2
    full.data <- rbind(old.data, new.data)
    
    #writing the data file
    utils::write.table(full.data,paste(filename,'_data.dat',sep=''),
      col.names=FALSE,row.names=FALSE,na='-9999',
      sep='\t',dec='.')
    
    parlabs <- unlist(regmatches(model, gregexpr('(?=\\().*?(?<=\\))', model, perl=T)))

    all.data <- rbind(new.data, old.data)
    
    results <- list(configural = NA, weak = NA, strong = NA, strict = NA)
    models <- list(configural = NA, weak = NA, strong = NA, strict = NA)
    
    for (invariance in names(results)) {
      equality <- character()
      if (invariance%in%c('weak','strong','strict')) equality <- 'lam'
      if (invariance%in%c('strong','strict')) equality <- paste(equality, 'alp', sep = '|')
      if (invariance%in%c('strict')) equality <- paste(equality, 'eps', sep = '|')
      
      parlabs_a <- parlabs
      if (invariance == 'configural') {
        parlabs_a <- gsub('\\)', 'A\\)', parlabs_a)
      } else {
        parlabs_a[!grepl(equality, parlabs_a)] <- gsub('\\)', 'A\\)', parlabs_a[!grepl(equality, parlabs_a)])
      }

      model_a <- paste(model, collapse = '\n')

      for (i in seq_along(parlabs)) {
        tmp_labs <- gsub('\\(', '\\\\(', parlabs)
        tmp_labs <- gsub('\\)', '\\\\)', tmp_labs)
        model_a <- gsub(tmp_labs[i], parlabs_a[i], model_a)
      }
      
      model_b <- gsub('A\\)', 'B\\)', model_a)
      
      analysis.options$model <- paste(gsub("\\(.+\\)", "", model), collapse = '\n')
      analysis.options$model <- paste(analysis.options$model, 'Model 1:', model_a, 'Model 2:', model_b, sep = '\n')
   
      args <- list(data=full.data,selected.items=selection$subtests,
        grouping=grouping,auxi=full.data[,NULL],suppress.model=TRUE,
        output.model=FALSE,svalues=FALSE,factor.structure=selection$parameters$factor.structure,
        filename=filename,cores=NULL,
        analysis.options=analysis.options)

      results[[invariance]] <- do.call('run.Mplus',args)
      results[[invariance]] <- as.data.frame(fitness(selection$parameters$objective, results[[invariance]], 'Mplus'))
      
      args$output.model <- TRUE
      models[[invariance]] <- do.call('run.Mplus',args)
    }
  }

  results <- do.call('rbind', results)
  if(any(sapply(full.data[, unlist(selection$subtests)], is.factor))) {
    warning('Model comparisons for ordinal indicators using Mplus are not yet implemented.', call. = FALSE)
  } else {
    results$`Chisq diff` <- NA
    results$`Df diff` <- NA
    results$`Pr(>Chisq)` <- NA
    
    # Model comparisons
    for (i in seq_along(models)[-1]) {
      m0 <- models[[i]]$summaries
      m1 <- models[[i-1]]$summaries
      correction <- ifelse(is.null(m0$ChiSqM_ScalingCorrection), 1, 
        (m0$Parameters * m0$ChiSqM_ScalingCorrection - m1$Parameters*m1$ChiSqM_ScalingCorrection)/(m0$Parameters - m1$Parameters))
      results$`Chisq diff`[i] <- -2*(m0$LL - m1$LL)/correction
      results$`Df diff`[i] <- m1$Parameters - m0$Parameters
      results$`Pr(>Chisq)`[i] <- stats::pchisq(results$`Chisq diff`[i], results$`Df diff`[i], lower.tail = FALSE)
    }
  }
  
  output <- list(comparison = results, models = models)
  return(output)
  
} # end function