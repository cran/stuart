### Roxygen-izable Documentation ----
#' Generating heuristics for the use in STUART subtest construction
#' 
#' Creates uninformative heuristic matrices for the use in \code{\link{mmas}}. 
#' 
#' This function generates a list of matrices which can be used as heuristics for all STUART constructions. This is mainly intended to write the structure of the heuristic matrices to an object, change components in line with theoretically derived heuristics and feed them back into \code{\link{mmas}} via the \code{heuristics} argument. The generated heuristics will contain only 1s and 0s, making it no heuristic information. Selection probabilities can be altered by manipulating the contents of the object created by \code{heuristics}. Setting a value to 0 will result in prohibiting a certain choice to be made. Please note, that it will lead to unpredictable behavior if the diagonal elements of the matrices produced in the arcs parameterization are set to values other than 0.
#' 
#' @author Martin Schultze
#' 
#' @seealso \code{\link{mmas}}
#' 
#' @concept ACO subtests
#' 
### Inputs ----
#' @param data A data.frame containing all relevant data.
#' @param factor.structure  A list linking factors to items. The names of the list elements correspond to the factor names. Each list element must contain a character-vector of item names that are indicators of this factor.
#' @param capacity A list containing the number of items per subtest. This must be in the same order as the \code{factor.structure} provided. If a single number, it is applied to all subtests. If \code{NULL} all items are evenly distributed among the subtests.
# #' @param number.of.subtests  A vector containing the number of subtests per construct. This must be in the same order as the \code{factor.structure} provided. If a single number, it is applied to all constructs. The default is to construct 1 subtest per construct.
#' @param repeated.measures A list linking factors that are repeated measures of each other. Repeated factors must be in one element of the list - other sets of factors in other elements of the list. When this is \code{NULL} (the default) a cross-sectional model is estimated.
#' @param mtmm A list linking factors that are measurements of the same construct with different methods. Measurements of the same construct must be in one element of the list - other sets of methods in other elements of the list. When this is \code{NULL} (the default) a single method model is estimated.
#' @param grouping The name of the grouping variable. The grouping variable must be part of \code{data} provided and must be a numeric variable.
#' @param localization Which parameterization to use when depositing pheromones. Can be either 'nodes' (the default) for depositing pheromones on selected nodes or 'arcs' for depositing on selection arcs.
#' @param ... Other arguments normally provided to \code{\link{mmas}}, which will be ignored.
#' 
### Outputs ---- 
#' @return Returns a list of the same length as the \code{factor.structure} argument provided.
#'   
#'   
### Examples ----
#' @examples
#' 
#' # heuristics for node localization
#' data(fairplayer)
#' fs <- list(si = names(fairplayer)[83:92])
#' 
#' (heu <- heuristics(fairplayer, fs, 4))
#' 
#' # Define anchor-item
#' heu$si[1] <- 10000
#' heu
#' 
#' # heuristics for arc localization
#' data(fairplayer)
#' fs <- list(si = names(fairplayer)[83:92])
#' 
#' (heu <- heuristics(fairplayer, fs, 4, localization = 'arcs'))
#' 
#' # Define equal selection of odd and even items
#' heu$si[1:10,] <- c(rep(c(0, 1), 5), rep(c(1, 0), 5))
#' heu
#' 
#' @export


### Function definition ----
heuristics <-
function(
  data, factor.structure, capacity=NULL, #number.of.subtests=1,
  repeated.measures=NULL, mtmm=NULL, grouping=NULL,  
  localization='nodes', ...
) { #begin function

  #combine arguments
  args <- as.list(match.call())[-1]
  args <- c(args,formals()[!names(formals())%in%names(args)])
  args$number.of.substests <- 1
  if (inherits(data, 'stuartHoldout')) {
    args$data <- data$calibrate
  }

  #sanity check
  objective <- NULL
  software <- 'lavaan'
  grouping <- NULL
  comparisons <- NULL
  do.call('sanitycheck',mget(names(formals(sanitycheck))))
  
  
  prepared <- do.call('data.prep',args)
  prepared$alpha <- 1

  #combine arguments
  args <- c(prepared,args[!names(args)%in%names(prepared)])
  args <- args[names(args)%in%names(formals(init.pheromones))]

  #generate pheromone matrices
  heuristics <- do.call('init.pheromones',args)

  #scale matrices
  heuristics <- lapply(heuristics, function(x) x^1/1e+100)

  class(heuristics) <- 'stuartHeuristics'
  attr(heuristics,'localization') <- localization
  return(heuristics)

} #end function
