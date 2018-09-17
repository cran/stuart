### Roxygen-izable Documentation ----
#' Subtest construction using a simple genetic algorithm
#' 
#' Construct subtests from a given pool of items using a simple genetic algorithm. Allows for multiple constructs, occasions, and groups.
#' 
### Details ----
#' The pheromone function provided via \code{objective} is used to assess the quality of the solutions. These functions can contain any combination of the fit indices provided by the estimation software. When using Mplus these fit indices are 'rmsea', 'srmr', 'cfi', 'tli', 'chisq' (with 'df' and 'pvalue'), 'aic', 'bic', and 'abic'. With lavaan any fit index provided by \code{\link[lavaan]{inspect}} can be used. Additionally 'crel' provides an aggregate of composite reliabilites, 'rel' provides a vector or a list of reliability coefficients for the latent variables, 'con' provides an aggregate consistency estimate for MTMM analyses, and 'lvcor' provides a list of the latent variable correlation matrices. For more detailed objective functions 'lambda', 'theta', 'psi', and 'alpha' provide the model-implied matrices. Per default a pheromone function using 'crel', 'rmsea', and 'srmr' is used. Please be aware that the \code{objective} must be a function with the required fit indices as (correctly named) arguments.
#' 
#' The genetic algorithm implemented selects parents in a two-step procedure. First, a fitness proportionate selection is performed to select \code{inviduals} times \code{reproduction} viable parents. Then, the non-self-adaptive version of mating proposed by Galán, Mengshoel, and Pinter (2013) is used to perform mating. In contrast to the original article, the \code{mating.index} and \code{mating.size} are handled as proportions, not integers. Similarity-based mating is based on the Jaccard Similarity. Mutation is currently always handled as an exchange of the selection state between two items. This results in mutation selecting one item that was not selected prior to mutation and dropping one item selected prior to mutation. Convergence is checked via the variance of the global-best values on the objective function, as proposed by Bhandari, Murthy, and Pal (2012). To avoid false convergence in the early search, the lower of either 10\% of the generations or 10 generations must be completed, before convergence is checked. For generalizability over different functions provided to \code{objective}, these are scaled to the first global-best found.
#' 
#' @author Martin Schultze
#' 
#' @seealso \code{\link{bruteforce}}, \code{\link{mmas}}, \code{\link{randomsamples}}
#' 
#' @concept ACO subtests
#' 
#' @references Bhandari, D., Murthy, C.A., & Pal, S.K. (2012). Variance as a Stopping Criterion for Genetic Algorithms with Elitist Model. Fundamenta Informaticae, 120, 145-164. doi:10.3233/FI-2012-754
#' @references Galán, S.F., Mengshoel, O.J., & Pinter,  R. (2013). A novel mating approach for genetic algorithms. Evolutionary Computation, 21(2), 197-229. doi:10.1162/EVCO_a_00067
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
#' @param auxiliary The names of auxiliary variables in \code{data}. These can be used in additional modeling steps that may be provided in \code{analysis.options$model}.
#' @param use.order A logical indicating whether or not to take the selection order of the items into account. Defaults to \code{FALSE}.
#' @param software The name of the estimation software. Can currently be 'lavaan' (the default) or 'Mplus'. Each option requires the software to be installed.
#' @param cores The number of cores to be used in parallel processing. If \code{NULL} (the default) the result of \code{\link[parallel]{detectCores}} will be used. On Unix-y machines parallel processing is implemented via \code{\link[parallel]{mclapply}}, on Windows machines it is realized via \code{\link[parallel]{parLapply}}.
#' @param objective A function that converts the results of model estimation into a pheromone. See 'details' for... details.
#' @param ignore.errors A logical indicating whether or not to ignore estimation problems (such as non positive-definite latent covariance matrices). Defaults to \code{FALSE}.
#' @param generations Maximum number of generations to run. Defaults to 128.
#' @param individuals The number of individuals per generation. Defaults to 64.
#' @param elitism The proportion of individuals from the last generation to carry over to the next generation. Defaults to 1/individuals, meaning that the best individual is retained into the next generation.
#' @param reproduction The proportion of individuals that are allowed to sire offspring. These individuals are selected using fitness proportionate selection. Defaults to .5. 
#' @param mutation The mutation probability. Defaults to .1. See 'details'.
#' @param mating.index The relative rank of the selected mate within the mating pool. A number bewteen 0, indicating a best-last mating, and 1 (the default), indicating a best-first mating. See 'details'.
#' @param mating.size The proportion of potential mates sampled from the pool of reproducers for each selected individual. Defaults to .25. See 'details'.
#' @param mating.criterion The criterion by which to select mates. Can be either 'similarity' or 'fitness' (the default). See 'details'.
#' @param tolerance The tolerance for deteriming convergence. Defaults to .0001. See 'details'.
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
#' \item{parameters}{A list of the parameters used.}
#' \item{analysis.options}{A list of the additional arguments passed to the estimation software.}
#' \item{timer}{An object of the class \code{proc_time} which contains the time used for the analysis.}
#' \item{log}{A \code{data.frame} containing the optimization history.}
#' \item{solution}{A list of matrices with the choices made in the global-best solution.}
#' \item{pheromones}{\code{NULL}}
#' \item{subtests}{A list containing the names of the selected items and their respective subtests.}
#' \item{final}{The results of the estimation of the global-best solution.}
#' 
#' 
### Examples ----
#' @examples
#' 
#' # Genetic selection in a simple situation
#' # requires lavaan
#' # number of cores set to 1 in all examples
#' data(fairplayer)
#' fs <- list(si = names(fairplayer)[83:92])
#' 
#' # minimal example
#' sel <- gene(fairplayer, fs, 4, 
#'   generations = 1, individuals = 10,  # minimal runtime, remove for application
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
#' # run to convergence
#' # switching to best-last mating and 50\% mating size
#' sel <- gene(fairplayer, fs, 4, 
#'   repeated.measures = repe, long.invariance = 'strong',
#'   mating.index = 0, mating.size = .5,
#'   seed = 55635, cores = 1)
#' 
#' # forcing a run through all generations
#' # by disabling the variance convergence rule
#' sel <- gene(fairplayer, fs, 4,
#'   repeated.measures = repe, long.invariance = 'strong',
#'   tolerance = 0, seed = 55635,
#'   cores = 1)
#' }
#' 
#' @export


### Function definition ----
gene <-
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
    
    auxiliary=NULL, use.order=FALSE,
    
    software='lavaan', cores=NULL,                                        #run settings
    
    objective=objective.preset, ignore.errors=FALSE,                      #fitness specs
    
    generations = 128, individuals = 64,                                  #algorithm specs
    elitism = 1/individuals, reproduction = .5, mutation = .1,
    mating.index = 1, mating.size = .25, 
    mating.criterion = 'fitness',
    tolerance = .0001,
    
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
    localization <- 'nodes'
    do.call('sanitycheck',mget(names(formals(sanitycheck))))
    
    timer <- proc.time()
    
    #check for software
    args$cores <- software.check(software,cores)
    
    #data preparation
    prepared <- do.call('data.prep',args)
    
    #combine arguments
    args <- c(prepared,args[!names(args)%in%names(prepared)])
    
    #call the core function
    solution <- do.call('stuart.gene',args)
    
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
    output$pheromones <- NULL
    output$subtests <- solution$selected.items
    output$final <- final.model
    
    class(output) <- 'stuartOutput'
    return(output)
    
  }