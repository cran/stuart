### Roxygen-izable Documentation ----
#' Generate an empirical objective function for item selection.
#' 
#' Generate an empirical objective function from default or empirical values for use in an item selection using STUART.
#' 
### Inputs ----
#' @param criteria A vector of names of criteria included in the objective function. Defaults to \code{c('rmsea', 'srmr', 'crel')}.
#' @param add A vector of names of criteria not used in the objective function, but added in order to be included in the log of solutions.
#' @param x Either a vector of values or an object of class \code{stuartOutput} from which to determine values in the objective function. If \code{NULL} (the default) values are generated from criteria-specific presets.
#' @param n Number of values to use in function determining. Defaults to 50, meaning if \code{side = 'top'} the 50 largest values are used to determine discrimination and difficulty parameters for each criterion.
#' @param side Which side good values are located at. \code{'top'} means large values are good (e.g. Reliability), \code{'bottom'} means small values are good (e.g. RMSEA), and \code{'middle'} means average values are good (e.g. factor correlations).
#' @param skew Whether to account for skew in the distribution using the \code{psn()} function from the \code{sn}-Package. Defaults to \code{FALSE}, meaning a normal distribution is used.
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
#' @seealso \code{\link{fixedobjective}}, \code{\link{extractobjective}}, \code{\link{objectivematrices}}
#' 
#' @export


empiricalobjective <- function(
  criteria = c('rmsea', 'srmr', 'crel'), 
  add = c('chisq', 'df', 'pvalue'),
  x = NULL,
  n = 50,
  side = NULL,
  skew = FALSE,
  scale = 1,
  matrices = NULL,
  fixed = NULL,
  comparisons = NULL,
  ...
) {

  called <- as.list(match.call())[-1]
  tmp <- as.list(formals(empiricalobjective))
  tmp <- tmp[setdiff(names(tmp), names(called))]
  called <- c(called, tmp)
  called <- called[names(called)!='x']
  
  if (is.null(x)) {
    warning('No empirical values provided to empiricalobjective. Attempting to return defaults.', call. = FALSE)
    out <- do.call(fixedobjective, called)
    class(out) <- 'stuartEmpiricalObjective'
    return(out)
  }
  
  if (inherits(x, 'stuartOutput')) {
    emp_values <- x$log[, criteria, drop = FALSE]
    emp_values <- as.list(emp_values)
  } else {
    tmp <- sapply(x, `[[`, 'solution.phe')
    tmp <- tmp[, !duplicated(t(tmp))]
    emp_values <- lapply(criteria, function(x) unlist(tmp[x, ]))
    names(emp_values) <- criteria
  }
  if (!is.null(matrices)) {
    cur_mat <- tmp[rownames(tmp) %in% names(matrices), , drop = FALSE]
    for (i in 1:nrow(cur_mat)) {
      emp_values[[length(emp_values) + 1]] <- cur_mat[i, ]
      names(emp_values)[length(emp_values)] <- rownames(cur_mat)[i]
    }
  }
  
  # predefined sets for typical criteria
  tops <- c('^(delta\\.)?([^c]+)rel', '^(delta\\.)?crel', '^(delta\\.)?cfi', '^(delta\\.)?tli', '^(delta\\.)?nnfi', '^(delta\\.)?rfi',
    '^(delta\\.)?nfi', '^(delta\\.)?pnfi', '^(delta\\.)?ifi', '^(delta\\.)?rni', '^(delta\\.)?logl', '^(delta\\.)?gfi', '^(delta\\.)?agfi', '^(delta\\.)?pgfi',
    '^(delta\\.)?mfi', '^(delta\\.)?ecvi', '^(delta\\.)?pvalue')
  tops <- paste0(tops, collapse = '|')
  bottoms <- c('^(delta\\.)?chisq', '^(delta\\.)?aic', '^(delta\\.)?bic', '^(delta\\.)?bic2', '^(delta\\.)?rmsea', '^(delta\\.)?rmr', '^(delta\\.)?srmr',
    '^(delta\\.)?crmr')
  bottoms <- paste0(bottoms, collapse = '|')
  mids <- c('^(delta\\.)?con')
  mids <- paste0(mids, collapse = '|')
  
  
  # check for and expand arguments
  end_reason <- NULL
  if (length(n) %in% c(1, length(criteria))) {
    n <- rep(n, length.out = length(criteria))
  } else {
    end_reason <- c(end_reason, 'n')
  }
  if (length(side) %in% c(1, length(criteria)) | is.null(side)) {
    if (is.null(side)) {
      side <- rep(NA, length(criteria))
      side[grep(tops, criteria)] <- 'top'
      side[grep(bottoms, criteria)] <- 'bottom'
      side[grep(mids, criteria)] <- 'center'
    }
    side <- rep(side, length.out = length(criteria))
  } else {
    end_reason <- c(end_reason, 'side')
  }
  if (length(skew) %in% c(1, length(criteria))) {
    skew <- rep(skew, length.out = length(criteria))
  } else {
    end_reason <- c(end_reason, 'skew')
  }
  if (length(scale) %in% c(1, length(criteria))) {
    scale <- rep(scale, length.out = length(criteria))
  } else {
    end_reason <- c(end_reason, 'scale')
  }
  if (!is.null(end_reason)) {
    addendum <- NULL
    if (any(end_reason=='scale') & !is.null(fixed)) {
      addendum <- ' Please scale components provided to \"fixed\" locally.'
    }
    stop(paste('Could not determine empirical objectives because arguments did not match the number of criteria. Problems with:', paste(end_reason, collapse = ', '), '.', addendum))
  }
  
  arguments <- data.frame(criteria, n, side, skew, scale)
  
  obj_list <- vector('list', length(criteria))
  names(obj_list) <- criteria
  for (i in criteria) {
    args <- as.list(subset(arguments, criteria == i))
    args$x <- emp_values[[i]]
    tmp <- do.call(extractobjective, args)
    tmp$string <- gsub('x', i, tmp$string)
    obj_list[[i]] <- tmp
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
          cur_values <- emp_values[[names(matrices)[i]]]
          cur_values <- cur_values[!is.na(cur_values)]
          args <- vector('list')
          if (length(js) > 1) {
            tmp <- lapply(cur_values, `[[`, j)
            args$side <- matrices[[i]]$side[[j]][k]
            args$skew <- matrices[[i]]$skew[[j]][k]
            args$scale <- matrices[[i]]$scale[[j]][k]
          } else {
            tmp <- cur_values
            args$side <- matrices[[i]]$side[k]
            args$skew <- matrices[[i]]$skew[k]
            args$scale <- matrices[[i]]$scale[k]
          }
          args$x <- sapply(tmp, `[`, k)
          args$n <- n[1]
          tmp <- do.call(extractobjective, args)
          if (length(js) > 1) {
            tmp$string <- gsub('x', paste0(names(matrices)[i], '[[', j, ']][', k, ']'), tmp$string)
          } else {
            tmp$string <- gsub('x', paste0(names(matrices)[i], '[', k, ']'), tmp$string)
          }
          obj_list[[length(obj_list) + 1]] <- tmp
        }
      }
    }
    
  }
  
  if (!is.null(fixed)) {
    if (inherits(fixed, 'function')) {
      fixed <- manualobjective(fixed)
    }
    obj_list[[length(obj_list) + 1]] <- fixed
    add <- union(add, eval(fixed$call$criteria))
    add <- union(add, eval(fixed$call$add))
  }
  
  
  tmp <- lapply(obj_list, `[[`, 'string')
  string <- paste0(tmp, collapse = ' + ')
  parsed <- parse(text = string)
  func <- function(...) eval(parsed)
  tmp <- vector('list', length(criteria) + length(add))
  names(tmp) <- c(criteria, add)
  formals(func) <- tmp
  
  out <- list(func = func, string = string, call = called)
  class(out) <- 'stuartEmpiricalObjective'
  return(out)
  
}
