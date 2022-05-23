### Roxygen-izable Documentation ----
#' Extracting empirical objective functions for item selection 
#' 
#' This is a bare-bones initial version of this approach.
#' 
### Inputs ----
#' @param x A vector of values for which to determine the objective function (e.g. RMSEA).
#' @param n  Number of values to use in function determining. Defaults to 50.
#' @param side Which side good values are located at. \code{'top'} means large values are good (e.g. Reliability), \code{'bottom'} means small values are good (e.g. RMSEA), and \code{'middle'} means average values are good (e.g. factor correlations).
#' @param skew Whether to account for skew in the distribution using the [sn::psn()] function. Defaults to \code{FALSE}, meaning a normal distribution is used.
#' @param scale A numeric scale to use in weighting the objective component. Defaults to 1.
#' @param ... Additional arguments.
#' 
#' @return Returns an object of class \code{stuartEmpiricalObjective}.
#' 
#' @author Martin Schultze
#' 
#' @seealso \code{\link{empiricalobjective}}, \code{\link{fixedobjective}}, \code{\link{objectivematrices}}
#' @keywords internal


### Empirical Objective Functions ----
extractobjective <- function(x,         # Input parameter
    n = 50,                             # Reference proportion
    side = c('top', 'bottom', 'center'),  # Where is good?
    skew = FALSE,                       # Use skew?
    scale = 1,                          # Scale output
    ...                                 # Additional arguments passed
  ) {
  
  # mark as deprecated
  .Deprecated('empiricalobjective', 'stuart', 'Objectives can now be extracted from stuartOutput-objects directly by using empiricalobjective().')
  
  # which side to use
  side <- side[1]
  if (side == 'centre') side <- 'center'
  
  # use only selected proportion
  y <- sort(x)
  ly <- length(y)
  py <- max(2, min(n, ceiling(length(y) * .25)))

  if (side == 'top') {
    y <- y[ly : (ly - py)]
  }
  if (side == 'bottom') {
    y <- y[1 : py]
  }
  if (side == 'center') {
    y <- y[((ly - py)/2) : ((ly + py) / 2)]
  }
  if (!(side %in% c('top', 'bottom', 'center'))) {
    stop('stop must be one of "top", "bottom", or "center".', .call = FALSE)
  }
  
  
  if (!skew) {
    m <- mean(y, na.rm = TRUE)
    s <- stats::sd(y, na.rm = TRUE)
    if (is.na(s) | s == 0) s <- .001
    if (side == 'top') {
      string <- paste0(scale, ' * pnorm(x, ', m, ', ', s, ', lower.tail = TRUE)')
    }
    if (side == 'bottom') {
      string <- paste0(scale, ' * pnorm(x, ', m, ', ', s, ', lower.tail = FALSE)')
    }
    if (side == 'center') {
      string <- paste0(scale, ' * 2 * ifelse(x > ', m, ', pnorm(x, ', m, ', ', s, ', lower.tail = FALSE), pnorm(x, ', m, ', ', s, ', lower.tail = TRUE))')
    }
  } else {
    sk_mod <- sn::selm(y ~ 1)
    pars <- sn::coef(sk_mod, 'DP')
    if (side == 'top') {
      string <- paste0(scale, ' * sn::psn(x, ', paste0(pars, collapse = ', '), ')')
    }
    if (side == 'bottom') {
      string <- paste0(scale, ' * (1 - sn::psn(x, ', paste0(pars, collapse = ', '), '))')
    }
    if (side == 'center') {
      string <- paste0(scale, ' * 2 * ifelse(x < ', sn::qsn(.5, dp = pars), ', sn::psn(x, ', paste0(pars, collapse = ', '), '), (1 - sn::psn(x, ', paste0(pars, collapse = ', '), ')))')
    }
    
  }
  
  parsed <- parse(text = string)
  func <- function(x) eval(parsed)
  
  out <- list(func = func, string = string)
  class(out) <- 'stuartEmpiricalObjective'
  return(out)
}