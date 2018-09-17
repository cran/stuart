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
#' @param invariance The invariance between the calibration and the validation sample. Can be one of 'configural', 'weak', 'strong', 'strict', or 'full', with the first being the default. Currently 'full' is only functional when using Mplus.
#' @param objective A function that converts the results of model estimation into a pheromone. If none is provided the default function \code{fitness} is used. This can be examined with \code{body(stuart:::fitness)}.
#' @param filename The stem of the filenames used to save inputs, outputs, and data files when \code{software='Mplus'}. This may include the file path. When \code{NULL} (the default) files will be saved to the temporary directory, which is deleted when the R session is ended.
#' 
### Outputs ----
#' @return Returns a list containing the \code{data.frame} \code{comparison} and an object containing the model results of the validation sample. 
#' 
#' \item{comparison}{A \code{data.frame} with 2 observations. The first observation shows the components of the objective function for the final model in the calibration sample. The second observation those of the model for the validation sample. Which variables are returned depends on the setting of \code{objective}.}
#' \item{validation}{When using \code{lavaan} for estimation, an object of class \code{lavaan} containing the model results fit to the validation sample. When using Mplus for estimation, a character vector containing the Mplus output for the validation sample.}
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
#' crossvalidate(sel, half1, half2, 
#'   invariance = 'strong')  # assuming equality of loadings and intercepts
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
#' crossvalidate(sel, split,
#'   invariance = 'weak')  # assuming equality of loadings
#' 
#' @export

### Function definition ----
crossvalidate <- 
function(
  selection, old.data, new.data,
  invariance = 'configural',
  objective = NULL,
  filename = NULL
) { #begin function
  
  if (!invariance%in%c('configural','weak','strong','strict','full'))
    stop('invariance must be configural, weak, strong, strict, or full.')
  
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

  # check fitness function
  if (is.null(objective)) objective <- selection$parameters$objective 
  args$objective <- objective

  # add analysis options
  args$analysis.options <- selection$analysis.options
  
  # run validation
  validated <- do.call(paste('crossvalidate',software,sep='.'),args)    
  fitness.options <- as.list(formals(fitness))
  fitness.options$solution.fit <- validated
  fitness.options$objective <- objective
  if ('con'%in%names(selection$log)) fitness.options$criteria <- c(as.character(fitness.options$criteria)[-1],'con')
  
  output <- list()
  output$comparison <- do.call(fitness,fitness.options)
  output$comparison <- rbind(selection$log[which.max(selection$log$pheromone),names(selection$log)%in%names(output$comparison)],array(data=unlist(output$comparison)))
  rownames(output$comparison) <- c('calibration','validation')
  
  args$output.model <- TRUE
  output$validation <- do.call(paste('crossvalidate',software,sep='.'),args)
  
  class(output) <- 'stuartCrossvalidate'
  return(output)
  
} #end function