### Roxygen-izable Documentation ----
#' Data selection for holdout validation.
#' 
#' Split a \code{data.frame} into two subsets for holdout validation.
#' 
#' @author Martin Schultze
#' 
#' @seealso \code{\link{crossvalidate}}
#' 
### Inputs ----
#' @param data A \code{data.frame}.
#' @param prop A single value or vector of proportions of data in calibration sample. Defaults to .5, for an even split.
#' @param grouping Name of the grouping variable. Providing a grouping variable ensures that the provided proportion is selected within each group.
#' @param seed A random seed. See \code{\link{Random}} for more details.
#' 
### Outputs ----
#' @return Returns a list containing two \code{data.frame}s, called calibrate and validate. The first corresponds to the calibration sample, the second to the validation sample.
#' 
#' @concept ACO subtests
#' 
### Examples ----
#' @examples
#' 
#' # seeded selection, 25% validation sample
#' data(fairplayer)
#' split <- holdout(fairplayer, .75, seed = 55635)
#' lapply(split, nrow) # check size of samples
#' 
#' @export


holdout <- function(data, prop = .5, grouping = NULL, seed = NULL) {
  
  # set random seed
  if (!is.null(seed)) {
    old.seed <- .Random.seed
    old.kind <- RNGkind()[1]
    set.seed(seed)
  }
  
  if (!is.null(grouping)) {
    if (length(prop) > 1 & length(prop) != nrow(unique(data[grouping]))) {
      prop <- prop[1]
      warning('The length of prop and the number of groups do not match. Only the first proportion is used.')
    }
    if (any(is.na(data[grouping]))) {
      data <- data[!is.na(data[grouping]), ]
      warning('Data contains observations with missing values on the grouping variables. These were excluded.')
    }
    n_cali <- ceiling(table(data[grouping]) * prop)
    filter <- NULL
    for (i in as.character(unlist(unique(data[grouping])))) {
      tmp <- which(data[grouping] == i)
      tmp_filter <- sample(tmp, n_cali[i])
      filter <- c(filter, tmp_filter)
    }
  } else {
    n_cali <- ceiling(nrow(data)*prop)
    filter <- sort(sample(nrow(data), n_cali))
  }
  
  output <- list(calibrate = data[filter, ], validate = data[-filter, ])
  class(output) <- 'stuartHoldout'
  
  # return to previous random seeds
  if (!is.null(seed)) {
    RNGkind(old.kind)
    .Random.seed <<- old.seed
  }
  
  return(output)
  
} 