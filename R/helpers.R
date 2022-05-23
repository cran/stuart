### Helper functions ----
# Fisher z
fishz <- function(x) log((1+x)/(1-x))/2

# inverse Fisher z
inv.fishz <- function(x) {
  out <- (exp(2*x)-1)/(exp(2*x)+1)
  out[is.na(out)] <- 1
  return(out)
}
