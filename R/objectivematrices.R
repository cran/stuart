### Roxygen-izable Documentation ----
#' Generate matrix-components for objective functions.
#' 
#' Generate objects of the correct structure for use in custom objective functions.
#' 
### Inputs ----
#' @param data A data.frame containing all relevant data.
#' @param factor.structure  A list linking factors to items. The names of the list elements correspond to the factor names. Each list element must contain a character-vector of item names that are indicators of this factor.
#' @param capacity A list containing the number of items per subtest. This must be in the same order as the \code{factor.structure} provided. If a single number, it is applied to all subtests. If \code{NULL} all items are evenly distributed among the subtests.
#' @param matrices Which matrix to extract. Can be one of \code{'lvcor'} (the default) for latent correlations, 'lambda', 'theta', 'psi', or 'alpha' for the model-implied matrices.
#' @param n.random The number of random draws to base values on. If 0 (the default) values are in the matrices are set to 0 and can be overwritten later. If any value larger than 0, the mean from \code{n.random} random solutions is used.
#' @param item.invariance A character vector of length 1 or the same length as \code{factor.structure} containing the desired invariance levels between items pertaining to the same subtest. Currently there are five options: 'congeneric', 'ess.equivalent', 'ess.parallel', 'equivalent', and 'parallel', the first being the default.
#' @param repeated.measures A list linking factors that are repeated measures of each other. Repeated factors must be in one element of the list - other sets of factors in other elements of the list. When this is \code{NULL} (the default) a cross-sectional model is estimated.
#' @param long.invariance A character vector of length 1 or the same length as \code{repeated.measures} containing the longitudinal invariance level of repeated items. Currently there are four options: 'configural', 'weak', 'strong', and 'strict'. Defaults to 'strict'. When \code{repeated.measures=NULL} this argument is ignored.
#' @param mtmm A list linking factors that are measurements of the same construct with different methods. Measurements of the same construct must be in one element of the list - other sets of methods in other elements of the list. When this is \code{NULL} (the default) a single method model is estimated.
#' @param mtmm.invariance A character vector of length 1 or the same length as \code{mtmm} containing the invariance level of MTMM items. Currently there are five options: 'none', 'configural', 'weak', 'strong', and 'strict'. Defaults to 'configural'. With 'none' differing items are allowed for different methods. When \code{mtmm=NULL} this argument is ignored.
#' @param grouping The name of the grouping variable. The grouping variable must be part of \code{data} provided and must be a numeric variable.
# #' @param group.invariance A single value describing the assumed invariance of subtests across groups. Currently there are four options: 'configural', 'weak', 'strong', and 'strict'. Defaults to 'strict'. When \code{grouping=NULL} this argument is ignored.
#' @param group.invariance A single value describing the assumed invariance of items across groups. Currently there are four options: 'configural', 'weak', 'strong', and 'strict'. Defaults to 'strict'. When \code{grouping=NULL} this argument is ignored.
#' @param comparisons A character vector containing any combination of 'item', 'long', 'mtmm', and 'group' indicating which invariance should be assessed via model comparisons. The order of the vector dictates the sequence in which model comparisons are performed. Defaults to \code{NULL} meaning that no model comparisons are performed. 
#' @param auxiliary The names of auxiliary variables in \code{data}. These can be used in additional modeling steps that may be provided in \code{analysis.options$model}.
#' @param use.order A logical indicating whether or not to take the selection order of the items into account. Defaults to \code{FALSE}.
#' @param software The name of the estimation software. Can currently be 'lavaan' (the default) or 'Mplus'. Each option requires the software to be installed.
#' @param cores The number of cores to be used in parallel processing. If \code{NULL} (the default) the result of \code{\link[parallel]{detectCores}} will be used. On Unix-y machines parallel processing is implemented via \code{\link[parallel]{mclapply}}, on Windows machines it is realized via \code{\link[parallel]{parLapply}}.
#' @param objective A function that converts the results of model estimation into a pheromone. See 'details' for... details.
#' @param ignore.errors A logical indicating whether or not to ignore estimation problems (such as non positive-definite latent covariance matrices). Defaults to \code{FALSE}.
#' @param analysis.options A list additional arguments to be passed to the estimation software. The names of list elements must correspond to the arguments changed in the respective estimation software. E.g. \code{analysis.options$model} can contain additional modeling commands - such as regressions on auxiliary variables.
#' @param suppress.model A logical indicating whether to suppress the default model generation. If \code{TRUE} a model must be provided in \code{analysis.options$model}.
#' @param ... Additional arguments passed either to \code{\link{randomsamples}} or to \code{lavaan}.
#' 
#' @return Returns an object of class \code{stuartFixedObjective}
#' 
#' @author Martin Schultze
#' 
#' @seealso \code{\link{empiricalobjective}}, \code{\link{extractobjective}}, \code{\link{objectivematrices}}
#' 
### Examples ----
#' @examples
#' 
#' # Extract latent correlation matric
#' # requires lavaan
#' # number of cores set to 1 in all examples
#' data(sups)
#' fs <- list(pro = names(sups)[2:13],
#'  fee = names(sups)[14:20])
#' 
#' mats <- objectivematrices(sups, fs, 3,
#'   cores = 1)
#' mats
#' 
#' mats$lvcor$use[,] <- FALSE
#' mats$lvcor$use[2, 1] <- TRUE
#' 
#' mats$lvcor$use
#' 
#' @export

objectivematrices <- 
  function(
    data, factor.structure, capacity=NULL, #number.of.subtests=1, #subtest settings
    matrices = c('lvcor'), n.random = 0,
    #invariance='parallel',
    item.invariance='congeneric',                  #cross invariance
    
    repeated.measures=NULL, long.invariance='strict', #long structure
    mtmm=NULL, mtmm.invariance='configural', #MTMM structure
    grouping=NULL, group.invariance='strict', #grouping structure
    comparisons=NULL,
    auxiliary=NULL, use.order=FALSE,
    software='lavaan', cores=NULL,                                        #run settings
    objective=NULL, ignore.errors=FALSE,                      #fitness specs
    analysis.options=NULL, suppress.model=FALSE,                          #modeling specs
    ...
  ) {
    
    if (is.null(objective)) {
      objective <- fixedobjective()
    }
    filt <-  matrices
    filt <- setdiff(filt,  c(eval(objective$call$criteria), eval(objective$call$add)))
    
    objective$call$add <- c(eval(objective$call$add), filt)
    
    addage <- vector('list', length(filt))
    names(addage) <- filt
    formals(objective$func) <- c(formals(objective$func), addage)
    
    args <- as.list(match.call())[-1]
    args <- c(args,formals()[!names(formals())%in%c(names(args),'...')])
    args <- args[names(args) %in% names(formals(randomsamples))]
    args$n <- max(c(1, n.random))
    args$objective <- objective
    args$matrices <- NULL
    
    attempts <- ifelse(n.random == 0, 0, 9)
    worked <- FALSE
    message('Attempting to extract matrices from a random subset...')
    while (!worked & attempts < 10) {
      invisible(utils::capture.output(suppressMessages(resi <- do.call(randomsamples, args))))
      attempts <- attempts + 1
      if (any(resi$log$pheromone != 0)) {
        worked <- TRUE
        message('\b done!')
      }
    }
    if (!worked) {
      stop('Was not able to extract matrices in 10 attempts. This may indicate a problem with the model, but may also be resolved by simply trying again.', call. = FALSE)
    }
    if (n.random == 0) {
      resi_mat <- resi$log_mat
    } else {
      tmp <- which(resi$log$pheromone != 0)[1]
      resi_mat <- lapply(resi$log_mat, `[`, tmp)
    }

    out <- vector('list', length(matrices))
    names(out) <- matrices
    types <- c('use', 'mean', 'sd', 'side', 'skew', 'scale')
    for (i in matrices) {
      out[[i]] <- vector('list', length(types))
      names(out[[i]]) <- types
      for (j in types) {
        out[[i]][[j]] <- extractmatrices(resi_mat[[i]][[1]], j)
      }
    }
    
    if (n.random > 0) {
      for (i in matrices) {
        dims <- dim(resi$log_mat[[i]][which(resi$log$pheromone != 0)[1]][[1]])
        if (i == 'lvcor') {
          resi_means <- do.call(rbind, sapply(resi$log_mat[[i]], c)) |> 
            fishz() |> colMeans(na.rm = TRUE) |> inv.fishz()
        } else {
          resi_means <- do.call(rbind, sapply(resi$log_mat[[i]], c)) |> colMeans(na.rm = TRUE)
        }
        resi_sds <- do.call(rbind, sapply(resi$log_mat[[i]], c)) |> apply(2, stats::sd, na.rm = TRUE)
        resi_sds <- sapply(resi_sds, \(x) max(c(.01, x)))

        out[[i]]$mean <- matrix(resi_means, nrow = dims[1], ncol = dims[2])
        out[[i]]$sd <-  matrix(resi_sds, nrow = dims[1], ncol = dims[2])
      }
    } 
    
    class(out) <- 'stuartObjectiveMatrices'
    return(out)

  }


extractmatrices <- function(x, type = 'use') {
  if (is.list(x)) {
    lapply(x, extractmatrices, type = type)
  }
  else {
    x <- unclass(x)
    if (type == 'use') {out <- x != 0}
    if (type == 'mean') {x[x != Inf] <- 0; out <- x}
    if (type == 'sd') {x[x != Inf] <- 1; out <- x}
    if (type == 'side') {x[x != Inf] <- 'center'; out <- x}
    if (type == 'skew') {out <- x == Inf}
    if (type == 'scale') {x[x != Inf] <- 1; out <- x}
    return(out)
  }
}
