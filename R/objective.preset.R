objective.preset <- function(chisq, df, pvalue, rmsea, srmr, crel){
    1 / (1 + exp(-10*(crel-.6))) +
    .5 * (1 - (1 / (1 + exp(-100*(rmsea-.05))))) +
    .5 * (1 - (1 / (1 +  exp(-100*(srmr-.06)))))
}

objective.preset.comparisons <- function(chisq, df, pvalue, rmsea, srmr, cfi, crel,
  delta.chisq, delta.df, delta.pvalue, delta.rmsea, delta.srmr, delta.cfi, delta.crel) {
  
  1 / (1 + exp(-10*(crel-.6))) +
  .5 * (1 - (1 / (1 + exp(-100*(rmsea-.05))))) +
  .5 * (1 - (1 / (1 +  exp(-100*(srmr-.06))))) +
  1 / (1 + exp(-30*(delta.crel-.1))) +
  .5 * (1 - (1 / (1 + exp(-300*(delta.rmsea-.01))))) +
  .5 * (1 - (1 / (1 + exp(-300*(delta.srmr-.01)))))    
}
