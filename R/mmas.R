### Roxygen-izable Documentation ----
#' Subtest construction using the Max-Min-Ant-System
#' 
#' Construct subtests from a given pool of items using the classical Max-Min Ant-System (Stützle, 1998). Allows for multiple constructs, occasions, and groups.
#' 
### Details ----
#' The pheromone function provided via \code{objective} is used to assess the quality of the solutions. These functions can contain any combination of the fit indices provided by the estimation software. When using Mplus these fit indices are 'rmsea', 'srmr', 'cfi', 'tli', 'chisq' (with 'df' and 'pvalue'), 'aic', 'bic', and 'abic'. With lavaan any fit index provided by \code{\link[lavaan]{inspect}} can be used. Additionally 'crel' provides an aggregate of composite reliabilites, 'rel' provides a vector or a list of reliability coefficients for the latent variables, 'con' provides an aggregate consistency estimate for MTMM analyses, and 'lvcor' provides a list of the latent variable correlation matrices. For more detailed objective functions 'lambda', 'theta', 'psi', 'alpha', and 'nu' provide the model-implied matrices. Per default a pheromone function using 'crel', 'rmsea', and 'srmr' is used. Please be aware that the \code{objective} must be a function with the required fit indices as (correctly named) arguments.
#' 
#' Using model comparisons via the \code{comparisons} argument compares the target model to a model with one less degree of assumed invariance (e.g. if your target model contains strong invariance, the comparison model contain weak invariance). Adding comparisons will change the preset for the objective function to include model differences. With comparisons, a custom objective function (the recommended approach) can also include all model fit indices with a preceding \code{delta.} to indicate the difference in this index between the two models. If more than one type of comparison is used, the argument of the objective function should end in the type of comparison requested (e.g. \code{delta.cfi.group} to use the difference in CFI between the model comparison of invariance across groups).
#' 
#' The scheduling of parameters is possible for the arguments \code{ants}, \code{colonies}, \code{evaporation}, \code{pbest}, \code{alpha}, \code{beta}, \code{tolerance}, and \code{deposit}. For all of these parameter scheduling is done when an array with two columns is provided. The first column of the array contains the timer, i.e. when to switch between parameter settings, the second column contains the values. The argument \code{schedule} can be used to select an absolute schedule (\code{schedule='run'}), a relative schedule which resets completely after a new global best is found (\code{schedule='colony'}), or a mixed version which resets the current phase of the schedule after a new global best is found (\code{schedule='mixed'}). When providing a parameter schedule for iterations 0, 3, and 10 using 'run' will result in a change after the third and the tenth iteration - irrespective of whether global best solutions were found. In contrast, using 'colony' will result in the first setting being used again once a new global best is found. This setting will then be used until iteration 3 (if no new best solution is found) before a switch occurs. If a new global best is found the setting will begin the sequence from the beginning. Using 'mixed' will result in the first setting being used until three consecutive iterations cannot produce a new global best. After this the second setting is used. If a new global best is found, the second setting is kept, but for the purpose of the schedule it is now iteration 3 again, meaning that the third setting will be used later than in a 'run' schedule.
#' 
#' @author Martin Schultze
#' 
#' @seealso \code{\link{bruteforce}}, \code{\link{gene}}, \code{\link{randomsamples}}, \code{\link{heuristics}}
#' 
#' @concept ACO subtests
#' 
#' @references Stützle, T. (1998). Local search algorithms for combinatorial problems: Analysis, improvements, and new applications. Unpublished doctoral dissertation. Darmstadt: Fachbereich Informatik, Universität Darmstadt.
#' 
### Inputs ----
#' @param data A data.frame containing all relevant data.
#' @param factor.structure  A list linking factors to items. The names of the list elements correspond to the factor names. Each list element must contain a character-vector of item names that are indicators of this factor.
#' @param capacity A list containing the number of items per subtest. This must be in the same order as the \code{factor.structure} provided. If a single number, it is applied to all subtests. If \code{NULL} all items are evenly distributed among the subtests.
#' @param item.weights A placeholder. Currently all weights are assumed to be one.
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
#' @param ants The number of ants per colony to be estimated. Can either be a single value or an array with two columns for parameter scheduling. See 'details'.
#' @param colonies The maximum number of colonies estimated since finding the latest global-best solution before aborting the process. Can either be a single value or an array with two columns for parameter scheduling. See 'details'.
#' @param evaporation The evaporation coefficient. Can either be a single value or an array with two columns for parameter scheduling. See 'details'.
#' @param alpha The nonlinearity coefficient of the pheromone-trail's contribution to determining selection probabilities. Defaults to 1 (linear). Can either be a single value or an array with two columns for parameter scheduling. See 'details'.
#' @param beta The nonlinearity coefficient of the heuristics' contribution to determining selection probabilities. Defaults to 1 (linear). Can either be a single value or an array with two columns for parameter scheduling. See 'details'.
#' @param pheromones A list of pheromones as created by \code{\link{mmas}}. This can be used to continue previous runs of this function.
#' @param heuristics An object of the class \code{stuartHeuristic} as provided by \code{\link{heuristics}} which contains heuristic information to be used in determining selection probabilities. If \code{NULL} (the default) selection probabilities are determined solely by the pheromones.
#' @param deposit Which deposit rule to use. Can be either 'ib' (the default) for an iteration-best deposit rule, or 'gb' for a global-best deposit rule.
#' @param localization Which localization to use when depositing pheromones. Can be either 'nodes' (the default) for depositing pheromones on selected nodes or 'arcs' for depositing on selection arcs.
#' @param pbest The desired overall probability of constructing the global-best solution when the algorithm convergels.  Can either be a single value or an array with two columns for parameter scheduling. See 'details'.
#' @param tolerance The tolerance of imprecision when comparing the pheromones to the upper and lower limits. Can either be a single value or an array with two columns for parameter scheduling. See 'details'.
#' @param schedule The counter which the scheduling of parameters pertains to. Can be either 'run' (the default), for a continuous schedule, 'colony', for a schedule that is restarted every time a new global best is found, or 'mixed' for a schedule that restarts its current phase every time a new global best is found. See 'details'.
#' @param analysis.options A list additional arguments to be passed to the estimation software. The names of list elements must correspond to the arguments changed in the respective estimation software. E.g. \code{analysis.options$model} can contain additional modeling commands - such as regressions on auxiliary variables.
#' @param suppress.model A logical indicating whether to suppress the default model generation. If \code{TRUE} a model must be provided in \code{analysis.options$model}.
#' @param seed A random seed for the generation of random samples. See \code{\link{Random}} for more details.
#' @param filename The stem of the filenames used to save inputs, outputs, and data files when \code{software='Mplus'}. This may include the file path. When \code{NULL} (the default) files will be saved to the temporary directory, which is deleted when the R session is ended.
#' 
#' 
#' 
### Outputs ---- 
#' @return Returns an object of the class \code{stuartOutput} for which specific \code{summary} and \code{plot} methods are available. The results are a list.
#' \item{call }{The called function.}
#' \item{software}{The software used to fit the CFA models.}
#' \item{parameters}{A list of the ACO parameters used.}
#' \item{analysis.options}{A list of the additional arguments passed to the estimation software.}
#' \item{timer}{An object of the class \code{proc_time} which contains the time used for the analysis.}
#' \item{log}{A \code{data.frame} containing the optimization history.}
#' \item{solution}{A list of matrices with the choices made in the global-best solution.}
#' \item{pheromones}{A list of matrices with the pheromones of each choice.}
#' \item{subtests}{A list containing the names of the selected items and their respective subtests.}
#' \item{final}{The results of the estimation of the global-best solution.}
#' 
### Examples ----
#' @examples
#' # MMAS in a simple situation
#' # requires lavaan
#' # number of cores set to 1 in all examples
#' data(fairplayer)
#' fs <- list(si = names(fairplayer)[83:92])
#' 
#' # minimal example
#' sel <- mmas(fairplayer, fs, 4, 
#'   colonies = 0, ants = 10,  # minimal runtime, remove for application
#'   seed = 55635, cores = 1)
#' summary(sel)
#' 
#' \donttest{
#' # longitudinal example
#' data(fairplayer)
#' fs <- list(si1 = names(fairplayer)[83:92],
#'   si2 = names(fairplayer)[93:102],
#'   si3 = names(fairplayer)[103:112])
#' 
#' repe <- list(si = c('si1', 'si2', 'si3'))
#' 
#' # change evaporation rate after 10 and 20 colonies
#' sel <- mmas(fairplayer, fs, 4, 
#'   repeated.measures = repe, long.invariance = 'strong',
#'   evaporation = cbind(c(0, 10, 20), c(.95, .8, .5)),
#'   seed = 55635, cores = 1)
#' 
#' # adding a predictor variable to selection model (using lavaan)
#' data(fairplayer)
#' fs <- list(si = names(fairplayer)[83:92])
#' 
#' add <- 'si ~ IGS'
#' 
#' sel <- mmas(fairplayer, fs, 4,
#'   auxiliary = 'IGS',
#'   analysis.options = list(model = add),
#'   seed = 55635, cores = 1)
#'
#' # inspect regression (in lavaan)
#' lavaan::coef(sel$final)
#' 
#' # same example, maximizing regression weight
#' obj <- function(chisq, df, pvalue, rmsea, srmr, beta) {
#'   beta[1, 'IGS']
#' }
#' 
#' sel <- mmas(fairplayer, fs, 4,
#'   auxiliary = 'IGS',
#'   analysis.options = list(model = add),
#'   objective = obj,
#'   seed = 55635, cores = 1)
#'   
#' # inspect regression (in lavaan)
#' lavaan::coef(sel$final)
#' }
#' 
#' @export


### Function definition ----
mmas <-
function(
  data, factor.structure, # data and structure
  
  capacity=NULL, item.weights=NULL,

  item.invariance='congeneric',                  #cross invariance

  repeated.measures=NULL,
  long.invariance='strict',
  
  mtmm=NULL,
  mtmm.invariance='configural',

  grouping=NULL,
  group.invariance='strict',

  comparisons=NULL,
  
  auxiliary=NULL, use.order=FALSE,

  software='lavaan', cores=NULL,                                        #run settings

  objective=NULL, ignore.errors=FALSE,                      #fitness specs

  ants=16, colonies=256, evaporation=.95,                               #general ACO specs
  alpha=1, beta=1, pheromones=NULL, heuristics=NULL,                    #general ACO specs
  deposit='ib', localization='nodes', pbest=.005, tolerance=.5,         #MMAS specs
  schedule='run',
  
  analysis.options=NULL, suppress.model=FALSE,                          #modeling specs
  seed=NULL,
  
  filename=NULL                                                         #stem of filenames for Mplus
) { #begin function

  #combine arguments
  args <- as.list(match.call())[-1]
  args <- c(args,formals()[!names(formals())%in%c(names(args),'...')])

  #select calibration sample (change to methods later)
  if (class(data) == 'stuartHoldout') {
    data <- data$calibrate
    args$data <- data
  }
  
  #sanity check
  do.call('sanitycheck',mget(names(formals(sanitycheck))))
  
  timer <- proc.time()

  #check for software
  args$cores <- software.check(software,cores)

  #data preparation
  prepared <- do.call('data.prep',args)

  #combine arguments
  args <- c(prepared,args[!names(args)%in%names(prepared)])

  #call the core function
  solution <- do.call('stuart.mmas',args)

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
  output$solution <- solution$solution.gb
  output$pheromones <- solution$pheromones
  output$subtests <- solution$selected.items
  output$final <- final.model$model

  class(output) <- 'stuartOutput'
  return(output)
  
}
