#' @export

plot.stuartOutput <-
function(x,remove.errors=TRUE,...) {
  
  phe <- x$log$pheromone
  run <- x$log$run
  
  if (remove.errors) {
    run <- run[phe!=0]
    phe <- phe[phe!=0]
  }
  
  best <- phe==max(phe)
  args <- as.list(match.call()[-1])
  args$x <- run
  args$y <- phe 
  args$xlab <- 'Run'
  args$ylab <- 'Pheromone'
  args$col <- as.numeric(best)+1
  if (is.null(args$pch)) args$pch <- 16

  do.call(graphics::plot,args)
  
  args$x <- stats::lowess(phe~run)
  args$y <- NULL
  do.call(graphics::lines,args)
  
  graphics::legend('bottomright',c('Alternative','Final Solution'),pch=args$pch,col=unique(args$col))
}
