### Roxygen-izable Documentation ----
#' Generate a fixed objective function for item selection.
#' 
#' Generate an objective function from default values for use in an item selection using STUART.
#' 
### Inputs ----
#' @param criteria A vector of names of criteria included in the objective function. Defaults to \code{c('rmsea', 'srmr', 'crel')}.
#' @param add A vector of names of criteria not used in the objective function, but added in order to be included in the log of solutions.
#' @param side Which side good values are located at. \code{'top'} means large values are good (e.g. Reliability), \code{'bottom'} means small values are good (e.g. RMSEA), and \code{'middle'} means average values are good (e.g. factor correlations).
#' @param scale A numeric scale to use in weighting the objective component. Defaults to 1.
#' @param matrices An object of class \code{stuartObjectiveMatrices} to include matrices (e.g. latent correlations) into the objective function.
#' @param fixed An object of class \code{stuartFixedObjective} to include already previously defined fixed objectives.
#' @param comparisons A vector of names of criteria included in the objective function which are related to model comparisons (e.g. when determining measurement invariance).
#' @param ... Additional arguments.
#' 
#' @return Returns an object of class \code{stuartFixedObjective}
#' 
#' @author Martin Schultze
#' 
#' @seealso \code{\link{empiricalobjective}}, \code{\link{extractobjective}}, \code{\link{objectivematrices}}
#' 
#' @export

fixedobjective <- function(
  criteria = c('rmsea', 'srmr', 'crel'), 
  add = c('chisq', 'df', 'pvalue'),
  side = NULL,
  scale = 1,
  matrices = NULL,
  fixed = NULL,
  comparisons = NULL,
  ...) {
  
  # predefined sets for typical criteria
  predef <- c('^([^c]+)rel', '^crel', '^cfi', '^tli', '^nnfi', '^rfi',
    '^nfi', '^pnfi', '^ifi', '^rni', '^gfi', '^agfi', '^pgfi',
    '^mfi', '^ecvi', '^pvalue', '^chisq', '^aic', '^bic', '^bic2', '^rmsea', '^rmr', '^srmr',
    '^crmr', '^con')
  predef_check <- paste0(predef, collapse = '|')
  
  if (!is.null(comparisons)) {
    predef_comp <- gsub('^\\^', '^delta.', predef)
    predef_comp[1] <- gsub('\\(\\[\\^c\\]\\+\\)', '', predef_comp[1])
    if (length(comparisons) > 1) {
      predef_comp <- outer(predef_comp, comparisons, paste, sep = '(.*).')
      predef_comp <- as.vector(predef_comp)
    }
    predef_check <- paste(predef_check, paste0(predef_comp, collapse = '|'),
      sep = '|')
  }
  
  if (any(!grepl(predef_check, criteria))) {
    end_reason <- criteria[!grepl(predef_check, criteria)]
    stop(paste0('Not all criteria provided to the objective function have defaults. Problem with: ', paste(end_reason, collapse = ', '), '. Consider using these as auxiliary information via \"add\" or including them manually via \"fixed\" instead.'))
  }
  
  defaults <- data.frame(criterion = predef, 
    side = c(rep('top', 16), rep('bottom', 8), rep('center', 1)),
    m = c(.7, .8, .95, .95, .95, .95, .95, .6, .95, .95, .95, .95, .5,
      .95, .4, .05, 0, 0, 0, 0, .05, .05, .05, .05, .7),
    s = c(.1, .075, .03, .03, .03, .03, .03, .1, .03, .03, .03, .03, .12,
      .03, .1, .1, 10, 10, 10, 10, .015, .02, .015, .02, .2),
    scale = 1)
  
  if (!is.null(comparisons)) {
    defaults_comp <- data.frame(criterion = predef_comp,
      side = defaults$side, 
      m = defaults$m / 10,
      s = defaults$s / 5,
      scale = defaults$scale)
    defaults_comp[grep('pvalue', defaults_comp$criterion), c('m', 's')] <- c(.05, .1) 
    defaults <- rbind(defaults, defaults_comp)
    if (!any(grepl('delta', c(criteria, add)))) {
      warning('Comparisons were requested, but the objective function does not contain any comparative indicators. Check whether this is intended - they can be added as \"criteria\" or via \"add\".', call. = FALSE)
    }
  }
  
  # Replace sides
  endreason <- vector('character')
  addendum <- vector('character')
  if (!is.null(side)) {
    if (length(side) %in% c(1, length(criteria))) {
      side <- rep(side, length.out = length(criteria))
    } else {
      endreason <- 'side'
    }
  }
  
  # Check for correct scaling
  if (length(scale) %in% c(1, length(criteria))) {
    scale <- rep(scale, length.out = length(criteria))
  } else {
    addendum <- NULL
    endreason <- 'scale'
    if (!is.null(fixed)) {
      addendum <- 'Please scale components provided to \"fixed\" beforehand.'
    }
  }

  if (length(endreason) > 0) {
    stop(paste0('Could not determine objectives because arguments did not match the number of criteria. Problems with: ', paste(endreason, sep = ', '), '. ', addendum))
  }
  
  if (is.null(side)) {
    tmp <- data.frame(criteria, scale)
    for (i in criteria) {
      filt <- sapply(defaults$criterion, grepl, x = i)
      defaults[filt, c('scale')] <- subset(tmp, criteria == i, select = c(scale))
    }
  } else {
    tmp <- data.frame(criteria, side, scale)
    for (i in criteria) {
      filt <- sapply(defaults$criterion, grepl, x = i)
      defaults[filt, c('side', 'scale')] <- subset(tmp, criteria == i, select = c(side, scale))
    }
  }
  
  obj_list <- vector('list', length = length(criteria))
  names(obj_list) <- criteria
  for (i in criteria) {
    filt <- sapply(defaults$criterion, grepl, x = i)
    cur_default <- defaults[filt, ]
    if (cur_default$side == 'top') {
      string <- paste0(cur_default$scale, ' * pnorm(x, ', cur_default$m, ', ', cur_default$s, ', lower.tail = TRUE)')
    }
    if (cur_default$side == 'bottom') {
      string <- paste0(cur_default$scale, ' * pnorm(x, ', cur_default$m, ', ', cur_default$s, ', lower.tail = FALSE)')
    }
    if (cur_default$side == 'center' | cur_default$side == 'centre') {
      string <- paste0(cur_default$scale, ' * 2 * ifelse(x > ', cur_default$m, ', pnorm(x, ', cur_default$m, ', ', cur_default$s, ', lower.tail = FALSE), pnorm(x, ', cur_default$m, ', ', cur_default$s, ', lower.tail = TRUE))')
    }
   
    string <- gsub('x', i, string)
    parsed <- parse(text = string)
    tmp_func <- function(x) eval(parsed)
    
    obj_list[[i]] <- list(func = tmp_func, string = string) 
  }
  
  if (!is.null(matrices)) {
    criteria <- c(criteria, names(matrices))
    for (i in 1:length(matrices)) {
      if (is.list(matrices[[i]]$use)) {
        js <- 1:length(matrices[[i]]$use)
      } else {
        js <- 1
      }
      for (j in js) {
        if (length(js) > 1) filt <- which(matrices[[i]]$use[[j]])
        else filt <- which(matrices[[i]]$use)
        for (k in filt) {
          if (length(js) > 1) {
            cur_side <- matrices[[i]]$side[[j]][k]
            cur_scale <- matrices[[i]]$scale[[j]][k]
            cur_m <-  matrices[[i]]$mean[[j]][k]
            cur_s <-  matrices[[i]]$sd[[j]][k]
          } else {
            cur_side <- matrices[[i]]$side[k]
            cur_scale <- matrices[[i]]$scale[k]
            cur_m <-  matrices[[i]]$mean[k]
            cur_s <-  matrices[[i]]$sd[k]
          }
          
          if (cur_side == 'top') {
            string <- paste0(cur_scale, ' * pnorm(x, ', cur_m, ', ', cur_s, ', lower.tail = TRUE)')
          }
          if (cur_side == 'bottom') {
            string <- paste0(cur_scale, ' * pnorm(x, ', cur_m, ', ', cur_s, ', lower.tail = FALSE)')
          }
          if (cur_side == 'center' | cur_side == 'centre') {
            string <- paste0(cur_scale, ' * 2 * ifelse(x > ', cur_m, ', pnorm(x, ', cur_m, ', ', cur_s, ', lower.tail = FALSE), pnorm(x, ', cur_m, ', ', cur_s, ', lower.tail = TRUE))')
          }
          if (length(js) > 1) {
            string <- gsub('x', paste0(names(matrices)[i], '[[', j, ']][', k, ']'), string)
          } else {
            string <- gsub('x', paste0(names(matrices)[i], '[', k, ']'), string)
          }
          parsed <- parse(text = string)
          tmp_func2 <- function(x) eval(parsed)
          
          tmp <- list(list(func = tmp_func2, string = string))
          names(tmp) <- paste0(names(matrices)[i], j, k)
          obj_list <- c(obj_list, tmp)
        }
      }
    }
    
  }
  
  if (!is.null(fixed)) {
    if (inherits(fixed, 'function')) {
      fixed <- manualobjective(fixed)
    }
    obj_list[[length(obj_list) + 1]] <- fixed
    add <- union(add, eval(names(formals(fixed$func))))
  }
  
  
  tmp <- lapply(obj_list, `[[`, 'string')
  string <- paste0(tmp, collapse = ' + ')
  parsed <- parse(text = string)
  func <- function(...) eval(parsed)
  tmp <- vector('list', length(criteria) + length(add))
  names(tmp) <- c(criteria, add)
  formals(func) <- tmp
  
  called <- as.list(match.call())[-1]
  tmp <- as.list(formals(empiricalobjective))
  tmp <- tmp[setdiff(names(tmp), names(called))]
  called <- c(called, tmp)
  called <- called[names(called)!='x']
  
  out <- list(func = func, string = string, call = called)
  class(out) <- 'stuartFixedObjective'
  return(out)
}


### Roxygen-izable Documentation ----
#' Convert empirical to fixed objective.
#' 
#' Convert an empirical objective to a fixed version to be used in item-selection. Sensible for extracting values from random selections and then using them in empirical but static objective functions.
#' 
### Inputs ----
#' @param x An object of class \code{stuartEmpiricalObjective}.
#' 
#' @return Returns an object of class \code{stuartFixedObjective}
#' 
#' @author Martin Schultze
#' 
#' @seealso \code{\link{empiricalobjective}}, \code{\link{fixedobjective}}
#' 
#' @export
as.stuartFixedObjective <- function(x) {
  if (!inherits(x, 'stuartEmpiricalObjective')) {
    warning('Only objects of class stuartEmpiricalObjective can be converted to stuartEmpiricalObjective. Input returned as is.')
  } else {
    class(x) <- 'stuartFixedObjective'
  }
  return(x)
}