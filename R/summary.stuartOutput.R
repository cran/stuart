#' @export

summary.stuartOutput <-
function(object,...) {
  Time <- as.numeric(object$timer[3])
  Models <- nrow(object$log)
  Replications <- sum(object$log$pheromone==max(object$log$pheromone))
  Results <- object$log[object$log$pheromone==cummax(object$log$pheromone),]
  Results <- Results[!duplicated(Results[,3:ncol(Results)]),]

  Out <- list(Subtests=object$subtests,Results=Results,Time=Time,Models=Models,
              Replications=Replications,Type=paste(object$call),Software=object$software)
  class(Out) <- 'summary.stuartOutput'
  return(Out)
}
