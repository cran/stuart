### Roxygen-izable Documentation ----
#' Subtest construction using a brute-force approach
#' 
#' Construct subtests from a given pool of items using a brute-force approach (i.e. by estimating all possible combinations).
#' 
### Details ----
#' The pheromone function provided via \code{objective} is used to assess the quality of the solutions. These functions can contain any combination of the fit indices provided by the estimation software. When using Mplus these fit indices are 'rmsea', 'srmr', 'cfi', 'tli', 'chisq' (with 'df' and 'pvalue'), 'aic', 'bic', and 'abic'. With lavaan any fit index provided by \code{\link[lavaan]{inspect}} can be used. Additionally 'crel' provides an aggregate of composite reliabilites, 'rel' provides a vector or a list of reliability coefficients for the latent variables, 'con' provides an aggregate consistency estimate for MTMM analyses, and 'lvcor' provides a list of the latent variable correlation matrices. For more detailed objective functions 'lambda', 'theta', 'psi', 'alpha', and 'nu' provide the model-implied matrices. Per default a pheromone function using 'crel', 'rmsea', and 'srmr' is used. Please be aware that the \code{objective} must be a function with the required fit indices as (correctly named) arguments.
#' 
#' Using model comparisons via the \code{comparisons} argument compares the target model to a model with one less degree of assumed invariance (e.g. if your target model contains strong invariance, the comparison model contain weak invariance). Adding comparisons will change the preset for the objective function to include model differences. With comparisons, a custom objective function (the recommended approach) can also include all model fit indices with a preceding \code{delta.} to indicate the difference in this index between the two models. If more than one type of comparison is used, the argument of the objective function should end in the type of comparison requested (e.g. \code{delta.cfi.group} to use the difference in CFI between the model comparison of invariance across groups).
#' 
#' @author Martin Schultze
#' 
#' @seealso \code{\link{mmas}}, \code{\link{gene}}, \code{\link{randomsamples}}, \code{\link{combinations}}
#' 
#' @concept ACO subtests
#' 
#' 
### Inputs ----
#' @param data A data.frame containing all relevant data.
#' @param factor.structure  A list linking factors to items. The names of the list elements correspond to the factor names. Each list element must contain a character-vector of item names that are indicators of this factor.
#' @param capacity A list containing the number of items per subtest. This must be in the same order as the \code{factor.structure} provided. If a single number, it is applied to all subtests. If \code{NULL} all items are evenly distributed among the subtests.
# #' @param number.of.subtests  A vector containing the number of subtests per construct. This must be in the same order as the \code{factor.structure} provided. If a single number, it is applied to all constructs. The default is to construct 1 subtest per construct.
#' @param item.invariance A character vector of length 1 or the same length as \code{factor.structure} containing the desired invariance levels between items pertaining to the same subtest. Currently there are five options: 'congeneric', 'ess.equivalent', 'ess.parallel', 'equivalent', and 'parallel', the first being the default.
#' @param repeated.measures A list linking factors that are repeated measures of each other. Repeated factors must be in one element of the list - other sets of factors in other elements of the list. When this is \code{NULL} (the default) a cross-sectional model is estimated.
#' @param long.invariance A character vector of length 1 or the same length as \code{repeated.measures} containing the longitudinal invariance level of repeated items. Currently there are four options: 'configural', 'weak', 'strong', and 'strict'. Defaults to 'strict'. When \code{repeated.measures=NULL} this argument is ignored.
#' @param mtmm A list linking factors that are measurements of the same construct with different methods. Measurements of the same construct must be in one element of the list - other sets of methods in other elements of the list. When this is \code{NULL} (the default) a single method model is estimated.
#' @param mtmm.invariance A character vector of length 1 or the same length as \code{mtmm} containing the invariance level of MTMM items. Currently there are five options: 'none', 'configural', 'weak', 'strong', and 'strict'. Defaults to 'configural'. With 'none' differing items are allowed for different methods. When \code{mtmm=NULL} this argument is ignored.
#' @param grouping The name of the grouping variable. The grouping variable must be part of \code{data} provided and must be a numeric variable.
#' @param group.invariance A single value describing the assumed invariance of items across groups. Currently there are four options: 'configural', 'weak', 'strong', and 'strict'. Defaults to 'strict'. When \code{grouping=NULL} this argument is ignored.
#' @param comparisons A character vector containing any combination of 'item', 'long', 'mtmm', and 'group' indicating which invariance should be assessed via model comparisons. The order of the vector dictates the sequence in which model comparisons are performed. Defaults to \code{NULL} meaning that no model comparisons are performed. 
#' @param auxiliary The names of auxiliary variables in \code{data}. These can be used in additional modeling steps that may be provided in \code{analysis.options$model}.
#' @param use.order A logical indicating whether or not to take the selection order of the items into account. Defaults to \code{FALSE}.
#' @param software The name of the estimation software. Can currently be 'lavaan' (the default), 'Mplus', or 'Mplus Demo'. Each option requires the software to be installed.
#' @param cores The number of cores to be used in parallel processing. If \code{NULL} (the default) the result of \code{\link[parallel]{detectCores}} will be used. On Unix-y machines parallel processing is implemented via \code{\link[parallel]{mclapply}}, on Windows machines it is realized via \code{\link[parallel]{parLapply}}.
#' @param objective A function that converts the results of model estimation into a pheromone. See \code{\link{mmas}} for details.
#' @param ignore.errors A logical indicating whether or not to ignore estimation problems (such as non positive-definite latent covariance matrices). Defaults to \code{FALSE}.
#' @param analysis.options A list additional arguments to be passed to the estimation software. The names of list elements must correspond to the arguments changed in the respective estimation software. E.g. \code{analysis.options$model} can contain additional modeling commands - such as regressions on auxiliary variables.
#' @param suppress.model A logical indicating whether to suppress the default model generation. If \code{TRUE} a model must be provided in \code{analysis.options$model}.
#' @param request.override The maximum number of combinations for which the estimation is performed immediately, without an additional override request.
#' @param filename The stem of the filenames used to save inputs, outputs, and data files when \code{software='Mplus'}. This may include the file path. When \code{NULL} (the default) files will be saved to the temporary directory, which is deleted when the R session is ended.
#' 
#' 
### Outputs ---- 
#' @return Returns an object of the class \code{stuartOutput} for which specific \code{summary} and \code{plot} methods are available. The results are a list.
#' \item{call }{The called function.}
#' \item{software}{The software used to fit the CFA models.}
#' \item{parameters}{A list of the ACO parameters used.}
#' \item{analysis.options}{A list of the additional arguments passed to the estimation software.}
#' \item{timer}{An object of the class \code{proc_time} which contains the time used for the analysis.}
#' \item{log}{A \code{data.frame} containing the estimation history.}
#' \item{solution}{\code{NULL}}
#' \item{pheromones}{\code{NULL}}
#' \item{subtests}{A list containing the names of the selected items and their respective subtests.}
#' \item{final}{The results of the estimation of the global-best solution.}
#' 
### Examples ----
#' @examples
#' 
#' # Bruteforce selection in a minimal example
#' # selecting 3 of 5 items
#' # requires lavaan
#' data(fairplayer)
#' fs <- list(ra = names(fairplayer)[53:57])
#' sel <- bruteforce(fairplayer, fs, 3,
#'   cores = 1)  # number of cores set to 1
#' summary(sel)  # Fit is perfect because of just-identified model
#'
#' @export


### Function definition ----
bruteforce <-
function(
  data, factor.structure, capacity=NULL, #number.of.subtests=1, #subtest settings

  #invariance='parallel',
  item.invariance='congeneric',                  #cross invariance

  repeated.measures=NULL, long.invariance='strict', #long structure

  mtmm=NULL, mtmm.invariance='configural', #MTMM structure
  
  grouping=NULL, group.invariance='strict', #grouping structure

  comparisons = NULL, auxiliary=NULL, use.order=FALSE,

  software='lavaan', cores=NULL,                                        #run settings

  objective=NULL, ignore.errors=FALSE,                      #fitness specs
  
  analysis.options=NULL, suppress.model=FALSE,                          #modeling specs

  request.override=10000,
  
  filename=NULL
) {#function begin

  #combine arguments
  args <- as.list(match.call())[-1]
  args <- c(args,formals()[!names(formals())%in%c(names(args),'...')])
  
  #select calibration sample (change to methods later)
  if (class(data) == 'stuartHoldout') {
    data <- data$calibrate
    args$data <- data
  }

  #sanity check
  localization <- 'nodes'
  do.call('sanitycheck',mget(names(formals(sanitycheck))))
  
  timer <- proc.time()

  #check for software
  args$cores <- software.check(software,cores)

  #data preparation
  prepared <- do.call('data.prep',args)

  #Feedbacking number of combinations
  combs <- do.call('compute.combinations', c(prepared[names(prepared)%in%names(formals(compute.combinations))],use.order))

  message(paste('There are',combs,'combinations that need to be tested.'))

  #ask for override if the estimation is not feasible
  if (combs>request.override) {
    override <- NA
    while(is.na(override)) {
      cat('\nWarning: Due to the number of possible combinations estimation may take very long or not be possible at all. Do you wish to continue?\n
 Enter y to continue or n to abort.\t'
      )
      override <- scan('',what='character',nmax=1,quiet=TRUE)
    }
    
    if (override!='y' & override!='Y') {
      stop('You aborted the estimation process.\n',call.=FALSE)
    }
  }

  #combine arguments
  args <- c(prepared,args[!names(args)%in%names(prepared)])

  solution <- do.call('stuart.bruteforce',args)

  args$data <- data
  args$output.model <- TRUE
  args$selected.items <- solution$selected.items
  args$selected <- solution$selected.gb

  tmp <- formals(paste('run',software,sep='.'))
  args <- args[names(args)%in%names(tmp)]
  args <- c(args,tmp[!names(tmp)%in%c(names(args))])

  final.model <- do.call(paste('run',software,sep='.'),args)

  #generating output
  output <- list(call=match.call())  
  output$software <- software
  output$parameters <- c(solution$parameters)
  output$analysis.options <- analysis.options
  output$timer <- proc.time() - timer
  output$log <- solution$log
  output$solution <- NULL
  output$pheromones <- NULL
  output$subtests <- solution$selected.items
  output$final <- final.model$model

  class(output) <- 'stuartOutput'
  return(output)

}
