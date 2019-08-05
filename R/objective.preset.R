objective.preset <- function(chisq, df, pvalue, rmsea, srmr, crel){
  1 / (1 + exp(6 - 10 * (crel))) +
    .5 * (1 - (1 / (1 + exp(5 - 100 * rmsea)))) +
    .5 * (1 - (1 / (1 +  exp(6 - 100 * srmr))))
}