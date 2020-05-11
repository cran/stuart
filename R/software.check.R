software.check <-
function(software,cores,...) { #function begin

  #check for Mplus if requested
  if (software=='Mplus') {
    if (!(nzchar(Sys.which('mplus'))|nzchar(Sys.which('/Applications/Mplus/mplus')))) {
      stop('Mplus is not installed. You could try using lavaan instead.\n',call.=FALSE) }
  }
  
  #check for lavaan if requested
  if (software=='lavaan') {
    if (!requireNamespace('lavaan')) {
      stop('lavaan is not installed. You could try or Mplus instead.\n',call.=FALSE) }
  }

  if (!requireNamespace('parallel')) {
    cat('parallel could not be loaded. STUART will continue without it.')
    cores <- 1
  }
  else {
    if (is.null(cores)) {
      cores <- parallel::detectCores()
    }
    if (cores > parallel::detectCores()) {
      cores <- parallel::detectCores()
    }    
  }  

  return(cores)

}
