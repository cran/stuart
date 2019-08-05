### Function definition ----
crossvalidate.lavaan <- 
function(
  selection, old.data, new.data, output.model=TRUE, 
  analysis.options = NULL, ...
) { #begin function
  
  # retrieve old results
  parameters <- lavaan::parTable(selection$final)
  parameters$est <- NULL
  parameters$se <- NULL
  parameters$start <- NULL

  if (!is.null(selection$call$grouping)) {
    stop('Crossvalidate is currently not available for situations with multiple groups.', call. = FALSE)
    # if (!selection$call$grouping %in% names(new.data) | 
    #     !all(unique(new.data[,selection$call$grouping])%in%unique(old.data[,selection$call$grouping]))) {
    #   warning('The validation sample contained a value on the grouping variable not contained in the calibration sample.',
    #     'Parameters of this group cannot be tested for invariance.')
    #   parameters <- parameters[parameters$group%in%c(0,1),]
    #   grouping <- NULL
    # } else {
    #   shrd <- which(unique(stats::na.omit(new.data[,selection$call$grouping]))%in%unique(stats::na.omit(old.data[,selection$call$grouping])))
    #   parameters <- parameters[parameters$group%in%c(shrd,0),]
    #   grouping <- selection$call$grouping
    #   if (length(unique(new.data[,grouping]))==1) {
    #     grouping <- NULL
    #     parameters$group <- NULL
    #   }
    # }
  } else {
    grouping  <- 'stuart_sample'
    old.data[, grouping] <- 'A'
    new.data[, grouping] <- 'B'
    
    par_1 <- parameters
    par_2 <- parameters
    par_2$group <- par_2$block <- 2
    n_par <- max(as.numeric(gsub('[^0-9]', '', par_1$plabel)), na.rm = TRUE)
    par_2$plabel[par_2$plabel != ''] <- paste0('.p', (1:n_par)+n_par, '.')
    par_2$free[par_2$free != 0] <- par_2$free[par_2$free != 0] + max(par_1$free)
    parameters <- rbind(par_1[par_1$plabel != '', ], par_2[par_2$plabel != '', ])
  }
  
  all.data <- rbind(new.data, old.data)

  results <- list(configural = NA, weak = NA, strong = NA, strict = NA)
  models <- list(configural = NA, weak = NA, strong = NA, strict = NA)
  
  for (invariance in names(results)) {
    pars <- parameters
    
    equality <- character()
    if (invariance%in%c('weak','strong','strict')) equality <- c(equality,'=~')
    if (invariance%in%c('strong','strict')) equality <- c(equality,'~1')
    if (invariance%in%c('strict')) equality <- c(equality,'~~')
    
    tmp <- pars[pars$label != '' & !pars$op%in%equality, ]
    pars[pars$label != '' & !pars$op%in%equality, 'label'] <- paste0(tmp$label, LETTERS[tmp$group])
    
    tmp <- pars[pars$label != '' & !duplicated(pars$label) & pars$free != 0, ]
    equal_labels <- lapply(tmp$label, function(x) pars$plabel[pars$label == x])
    par_constraints <- pars[FALSE, ]
    if (length(unlist(lapply(equal_labels, function(x) x[-1]))) > 0) {
      par_constraints[1:length(unlist(lapply(equal_labels, function(x) x[-1]))), ] <- NA
      
      par_constraints$lhs <- unlist(lapply(equal_labels, function(x) rep(x[1], length(x)-1)))
      par_constraints$op <- '=='
      par_constraints$rhs <- unlist(lapply(equal_labels, function(x) x[-1]))
      par_constraints$user <- 2
      par_constraints[, c('block', 'group', 'free', 'exo')] <- 0
      par_constraints[, c('label', 'plabel')] <- ''
    }
    
    pars <- rbind(pars, par_constraints)
    pars$id <- 1:nrow(pars)
    
    analysis.options$model <- pars
    
    args <- list(data=all.data,selected.items=selection$subtests,
      grouping=grouping,auxi=all.data[,NULL],suppress.model=TRUE,
      analysis.options=analysis.options,objective=selection$parameters$objective,ignore.errors=TRUE,
      output.model=FALSE,factor.structure=selection$parameters$factor.structure)
    
    results[[invariance]] <- do.call('run.lavaan',args)
    results[[invariance]] <- as.data.frame(fitness(selection$parameters$objective, results[[invariance]], 'lavaan'))
    
    args$output.model <- TRUE
    models[[invariance]] <- do.call('run.lavaan',args)
    
    
  }

  results <- do.call(rbind, results)
  comps <- do.call(lavaan::lavTestLRT, c(object=models[[1]], models[-1]))
  rownames(comps) <- names(models)
  results <- cbind(results, comps[, 5:7])
  output <- list(comparison = results, models = models)
  
  return(output)
  
} #end function