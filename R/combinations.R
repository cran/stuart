### Roxygen-izable Documentation ----
#' Compute the number of possible subtest combinations
#' 
#' Used to compute the number of possible subtest constellations prior to performing item selection.
#' 
#' 
#' @author Martin Schultze
#' 
#' @seealso \code{\link{bruteforce}}, \code{\link{mmas}}, \code{\link{gene}}
#' 
#' @concept ACO subtests
#' 
### Inputs ----
#' @param data A data.frame containing all relevant data.
# #' @param number.of.subtests  A vector containing the number of subtests per construct. This must be in the same order as the \code{factor.structure} provided. If a single number, it is applied to all constructs. The default is to construct 1 subtest per construct.
#' @param factor.structure  A list linking factors to items. The names of the list elements correspond to the factor names. Each list element must contain a character-vector of item names that are indicators of this factor.
#' @param capacity A list containing the number of items per subtest. This must be in the same order as the \code{factor.structure} provided. If a single number, it is applied to all subtests. If \code{NULL} all items are evenly distributed among the subtests.
#' @param repeated.measures A list linking factors that are repeated measures of each other. Repeated factors must be in one element of the list - other sets of factors in other elements of the list. When this is \code{NULL} (the default) a cross-sectional model is estimated.
#' @param mtmm A list linking factors that are measurements of the same construct with different methods. Measurements of the same construct must be in one element of the list - other sets of methods in other elements of the list. When this is \code{NULL} (the default) a single method model is estimated.
#' @param use.order A logical indicating whether or not to take the selection order of the items into account. Defaults to \code{FALSE}.
#' @param ... Other arguments normally provided to \code{\link{mmas}}, which will be ignored.
#' 
### Outputs ---- 
#' @return Returns the number of possible subtest constellations.
#' 
### Examples ----
#' @examples
#' 
#' # Determine number of combinations in a simple situation
#' data(fairplayer)
#' fs <- list(si = names(fairplayer)[83:92])
#' combinations(fairplayer, fs, 4)
#' 
#' # Number of combinations with repeated measures
#' data(fairplayer)
#' fs <- list(si1 = names(fairplayer)[83:92],
#'   si2 = names(fairplayer)[93:102],
#'   si3 = names(fairplayer)[103:112])
#' repe <- list(si = c('si1', 'si2', 'si3'))
#' combinations(fairplayer, fs, 4, repeated.measures = repe)
#' 
#' @export
#' 


### Function definition ----
combinations <-
function(
  data, factor.structure, capacity=NULL, #subtest settings
  repeated.measures=NULL, mtmm=NULL, use.order=FALSE,...
) {#function begin

  #arguments
  args <- as.list(match.call())[-1]
  args <- c(args,formals()[!names(formals())%in%c(names(args),'...')])
  args <- c(args,formals(data.prep)[!names(formals(data.prep))%in%c(names(args),'...')])
  args$number.of.subtests <- 1
  
  #sanity check
  objective <- NULL
  localization <- 'nodes'
  do.call('sanitycheck',mget(names(formals(sanitycheck))))
  
  #data preparation
  if (class(data) == 'stuartHoldout') {
    args$data <- data$calibrate
  }
  prepared <- do.call('data.prep',args)

  combs <- do.call('compute.combinations', c(prepared[names(prepared)%in%names(formals(compute.combinations))],use.order))

  return(combs)
}
