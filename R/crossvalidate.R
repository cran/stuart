### Roxygen-izable Documentation ----
#' Cross-Validate a Measurement Model
#' 
#' Cross-validate a measurement model obtained from STUART.
#' 
#' @author Martin Schultze
#' 
#' @seealso \code{\link{holdout}}, \code{\link{mmas}}, \code{\link{bruteforce}}
#' 
### Inputs ----
#' @param selection An object of class \code{stuartOutput}.
#' @param old.data A \code{data.frame} of the calibration sample.
#' @param new.data A \code{data.frame} of the validation sample.
#' @param filename The stem of the filenames used to save inputs, outputs, and data files when \code{software='Mplus'}. This may include the file path. When \code{NULL} (the default) files will be saved to the temporary directory, which is deleted when the R session is ended.
#' 
### Outputs ----
#' @return Returns a list containing the \code{data.frame} \code{comparison} and an object containing the model results of the four different invariance assumptions. 
#' 
#' \item{comparison}{A \code{data.frame} with 4 observations, each observation representing a level of measurement invariance. The number of columns depends on the arguments of the \code{objective} used in the original selection. In addition to those columns, three additional columns with the (corrected) Likelihood-Ratio-Tests are reported.}
#' \item{models}{A list of the four model results either of class \code{lavaan} or \code{mplus.model}, depending on the \code{software}-setting of the original selection.}
#' 
#' @concept ACO subtests
#' 
### Examples ----
#' @examples
#' 
#' # Split data into two halves
#' data(fairplayer)
#' half1 <- fairplayer[1:72,]
#' half2 <- fairplayer[73:143,]
#' 
#' # Simple example from bruteforce
#' fs <- list(ra = names(fairplayer)[53:57])
#' sel <- bruteforce(half1, fs, 3,
#'   cores = 1)  # number of cores set to 1
#' 
#' # Validation
#' crossvalidate(sel, half1, half2)
#' 
#' # Using the 'holdout' function for data split
#' data(fairplayer)
#' split <- holdout(fairplayer, seed = 55635)
#' 
#' # Simple example from bruteforce
#' fs <- list(ra = names(fairplayer)[53:57])
#' sel <- bruteforce(split, fs, 3,
#'   cores = 1)  # number of cores set to 1
#' 
#' # Validation
#' crossvalidate(sel, split)
#' 
#' @export

### Function definition ----
crossvalidate <- 
function(
  selection, old.data, new.data,
  filename = NULL
) { #begin function
  
  # check estimation software
  software <- selection$software
  
  if (software=='Mplus' & is.null(old.data)) stop('When using Mplus the old.data is required.')
  
  # set arguments for software-specific runs
  args <- as.list(match.call())[-1]
  args <- c(args,formals()[!names(formals())%in%c(names(args),'...')])
  args <- args[names(args)%in%names(formals(paste('crossvalidate',software,sep='.')))]
  
  #select calibration sample (change to methods later)
  if (class(old.data) == 'stuartHoldout') {
    args$new.data <- old.data$validate
    args$old.data <- old.data$calibrate
  }

  # add analysis options
  args$analysis.options <- selection$analysis.options
  
  # run validation
  output <- do.call(paste('crossvalidate',software,sep='.'),args)
  
  class(output) <- c('stuartCrossvalidate')
  return(output)
  
} #end function