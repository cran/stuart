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
#' @param determined Name of a variable indicating the pre-determined assignment to the calibration or the validation sample. This variable must be a factor containing only \code{NA} (no determined assingment), \code{"calibrate"}, or \code{"validate"}. If no variable is provided (the default) all cases are assigned randomly.
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


holdout <- function(data, prop = .5, grouping = NULL, seed = NULL,
  determined = NULL) {
  
  # set random seed
  if (!is.null(seed)) {
    old.seed <- .Random.seed
    old.kind <- RNGkind()[1]
    set.seed(seed)
  }
  
  if (is.null(determined)) {
    determined <- 'determined'
    data[, determined] <- NA
  } else {
    if (any(!is.na(data[determined]) & !(unlist(data[determined]) %in% c('validate', 'calibrate')))) 
      stop('The variable provided to "determined" contained values other than NA, "validate" or "calibrate".', call. = FALSE)
  }

  data[, determined] <- factor(data[, determined], 
    levels = c('calibrate', 'validate'))
  det_cali <- which(data[determined] == 'calibrate')
  det_vali <- which(data[determined] == 'validate')

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
    tmp <- as.factor(data[, grouping])
    n_det <- table(tmp[data[determined] == 'calibrate'])
    n_cali <- n_cali - n_det

    if (any(n_cali < 0)) {
      n_cali[n_cali < 0] <- 0
      warning('The deterministic sample assignment exceeded the proportion provided to "prop". Resulting proportions may differ.')
    }
    
    filter <- NULL
    for (i in as.character(unlist(unique(data[grouping])))) {
      tmp <- which(data[grouping] == i)
      tmp <- tmp[!(tmp %in% c(det_cali, det_vali))]
      if (length(tmp) > 0) {
        tmp_filter <- sample(tmp, min(n_cali[i], length(tmp)))
      } else {
        tmp_filter <- 0
        warning('The deterministic assignment lead to unbalanced samples. Please check the results carefully.')
      }
      filter <- c(filter, tmp_filter)
    }
  } else {
    n_cali <- ceiling(nrow(data)*prop)
    n_cali <- n_cali - sum(data[determined] == 'calibrate', na.rm = TRUE)
    if (n_cali < 0) {
      n_cali <- 0
      warning('The deterministic sample assignment exceeded the proportion provided to "prop". Resulting proportions may differ.')
    }
    tmp <- 1:nrow(data)
    tmp <- tmp[!(tmp %in% c(det_cali, det_vali))]
    if (length(tmp) > 0) {
      filter <- sort(sample(tmp, min(n_cali, length(tmp))))
    } else {
      filter <- 0
      warning('The deterministic assignment lead to unbalanced samples. Please check the results carefully.')
    }
  }
  filter <- c(det_cali, filter)
  
  output <- list(calibrate = data[filter, ], validate = data[-filter, ])
  class(output) <- 'stuartHoldout'
  
  # return to previous random seeds
  if (!is.null(seed)) {
    RNGkind(old.kind)
    .Random.seed <<- old.seed
  }
  
  return(output)
  
} 