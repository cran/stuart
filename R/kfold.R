### Roxygen-izable Documentation ----
#' k-Folds Crossvalidation
#' 
#' k-Folds crossvalidation for item selection using any approach implemented in STUART.
#' 
#' @details The function splits the provided data into k subsets using \code{\link{holdout}} and runs the item-selection procedure requested via \code{type} on the training datasets. Validation is performed using \code{\link{crossvalidate}} to check for invariance of the measurement models between the training and validation data up to the invariance level provided via \code{max.invariance}. The final item selection is based on the highest value on the objective function in the multiple-group SEM imposing \code{max.invariance} between the training and validation data.
#' 
#' @author Martin Schultze
#' 
#' @seealso \code{\link{holdout}} \code{\link{crossvalidate}}
#' 
#' @param type A \code{character} calling the item-selection procedure. Can be one of "mmas" (see \code{\link{mmas}}), "gene" (see \code{\link{gene}}), "bruteforce" (see \code{\link{bruteforce}}), or randomsamples (see \code{\link{randomsamples}}).
#' @param k The number of folds.
#' @param max.invariance The maximum measurement invariance level which will be tested. Currently there are four options: 'configural', 'weak', 'strong', and 'strict' (the default). All levels below \code{max.invariance} are also tested. 
#' @param seed The random seed. 
#' @param seeded.search A \code{logical} indicating whether the \code{seed} should also be used for the search procedure (the default) or only for the sample splitting.
#' @param ... Arguments passed to the item-selection procedure called with \code{type}.
#' @param remove.details A \code{logical} indicating whether to remove detailed information such as models and copies of datasets. Reduces output size by approx. 90\%. Defaults to \code{TRUE}.
#' 
#' @return Returns an object of the class \code{stuartKfold} for which specific \code{summary} and \code{print} methods are available. The results are a list.
#' \item{call }{The called function.}
#' \item{subtests}{A list containing the names of the selected items and their respective subtests.}
#' \item{solution}{A list of matrices with the choices made in the global-best solution.}
#' \item{final}{The results of the estimation of the solution leading to best objective value when cross-validated.}
#' \item{frequencies}{A list of matrices showing the relative frequencies with which an item was selected across the k-folds.}
#' \item{full}{A list of the results returned by the k runs of \code{type}.}
#' \item{crossvalidations}{A list of data.frames showing the fit and model comparisons of all invariance levels up to \code{max.invariance} in each of the k folds.}
#' \item{data}{A \code{data.frame}. The same as the original \code{data.frame} provided to \code{data} with the added variable \code{stuartKfold} indicating which fold an observation was assigned to.}
#' 
#' @examples
#' 
#' # k-Folding for a simple bruteforce selection
#' data(fairplayer)
#' fs <- list(ra = names(fairplayer)[53:57])
#' 
#' sel <- kfold('bruteforce', k = 2,
#'   data = fairplayer, factor.structure = fs, 
#'   capacity = 3, seed = 55635,
#'   cores = 1)
#' summary(sel)
#' 
#' \donttest{
#' ### longitudinal example with mmas ----
#' data(fairplayer)
#' fs <- list(si1 = names(fairplayer)[83:92],
#'   si2 = names(fairplayer)[93:102],
#'   si3 = names(fairplayer)[103:112])
#' 
#' repe <- list(si = c('si1', 'si2', 'si3'))
#' 
#' sel_mmas <- kfold('mmas', k = 3,
#'   data = fairplayer, factor.structure = fs,
#'   repeated.measures = repe, long.invariance = 'strong',
#'   capacity = 3, seed = 55635, pbest = .5,
#'   cores = 1)
#' summary(sel_mmas)
#' }
#' 
#' 
#' @export


kfold <- function(type, k = 5,
  max.invariance = 'strict',
  seed = NULL, seeded.search = TRUE,
  ...,
  remove.details = TRUE) {

  # unpack ellipses
  args <- list(..., seed = seed)
  
  # check for multiple groups
  if ('grouping' %in% names(args)) stop('Multiple groups are not yet supported in k-folds crossvalidation.', call. = FALSE)
  
  # split data
  folded <- list()
  
  hold_args <- args[names(args) %in% names(formals(holdout))]
  hold_args$data <- args$data
  hold_args$data$.determined_internal <- NA
  hold_args$prop <- 1/k
  hold_args$determined <- '.determined_internal'
  
  for (i in 1:k) {
    folded[[i]] <- do.call('holdout', hold_args)
    hold_args$data <- do.call('rbind', folded[[i]])
    hold_args$data$.determined_internal[1:nrow(folded[[i]]$calibrate)] <- 'validate'
    tmp <- folded[[i]]
    folded[[i]]$calibrate <- tmp$validate
    folded[[i]]$validate <- tmp$calibrate
  }
  
  # Run searches
  searches <- list()

  run_args <- args[names(args) %in% names(formals(type))]
  
  if (!seeded.search) run_args$seed <- NULL
  
  for (i in 1:k) {
    run_args$data <- folded[[i]]
    message(paste0('\nRunning fold number ', i, ' of ', k, '.\n'))
    searches[[i]] <- try(do.call(type, run_args))
    if ('try-error' %in% class(searches[[i]])) {
      warning(paste0('The search procedure did not terminate normally in fold ', i), call. = FALSE)
    }
  }
  
  check <- sapply(searches, function(x) 'try-error' %in% class(x))
  if (all(check)) stop('None of the folds resulted in viable solutions. This may be the result of the sample being too small for the number of folds.', call. = FALSE)

  # Run crossvalidation
  message('\nRunning cross-validation.\n')
  cv <- list()
  for (i in 1:k) {
    selection <- searches[[i]]
    old.data <- folded[[i]]
    # invisible(capture.output(cv[[i]] <- suppressWarnings(try(crossvalidate(selection, old.data, max.invariance = max.invariance), silent = TRUE))))
    cv[[i]] <- suppressWarnings(try(crossvalidate(selection, old.data, max.invariance = max.invariance)))
    if ('try-error' %in% class(cv[[i]])) {
      cv[[i]] <- list(comparison = NULL, models = NULL, matrices = NULL)
      warning(paste0('The crossvalidation produced an error in fold ', i), call. = FALSE)
    }
  }
  
  # Reorganize solutions
  tmp <- sapply(searches, function(x) 'try-error' %in% class(x))
  solu <- searches[[which.min(tmp)]][['solution']]
  for (i in seq_along(solu)) {
    solu[[i]] <- do.call('rbind', lapply(lapply(searches[!tmp], `[[`, 'solution'), `[[`, i))
  }
  
  # identify best solution
  phe <- sapply(cv, function(x) {
    if (is.null(x[['comparison']])) return(NA)
    else return(x[['comparison']][max.invariance, 'pheromone'])})
  
  best <- searches[[which.max(phe)]]
  best_cv <- cv[[which.max(phe)]]
  
  # remove models and data for non-best
  if (remove.details) {
    searches <- lapply(searches, function(x) x[!(names(x) %in% c('final', 'call'))])
    cv <- lapply(cv, `[[`, 'comparison')
  }
  
  dats <- do.call(rbind, lapply(folded, `[[`, 'validate'))
  dats$stuartKfold <- unlist(lapply(1:k, function(x) rep(x, nrow(folded[[x]][['validate']]))))
  rownames(dats) <- NULL
  
  out <- list(call = match.call())
  out$subtests <- best$subtests
  out$solution <- best$solution
  out$final <- best_cv[['models']][[max.invariance]]
  out$validation <- best_cv[['comparison']]
  out$frequencies <- lapply(solu, colMeans)
  out$full <- searches
  out$crossvalidations <- cv
  out$data <- dats
  
  class(out) <- 'stuartKfold'
  
  return(out)

}