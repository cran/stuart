compute.comparisons <- function(objective, baseline, alternative, type) {
  
  deltas <- grep('delta\\.', names(formals(objective)), value = TRUE)
  tmp_filt <- sapply(strsplit(deltas, '\\.'), length)==2
  deltas[tmp_filt] <- paste0(deltas, '.', type)
  deltas <- grep(paste0('\\.',type), deltas, value = TRUE)
  
  tmp_deltas <- sapply(strsplit(deltas, '\\.'), function(x) x[2])

  if (length(tmp_deltas)==0 | length(baseline) == 1 | length(alternative) == 1) {
    return(tmp_deltas)
  }
  
  diffs <- as.list(unlist(alternative[tmp_deltas])-unlist(baseline[tmp_deltas]))
  names(diffs) <- deltas
  
  if (any(grepl('pvalue|chisq', tmp_deltas))) {

    if (class(baseline$model)[1]=='lavaan') {
      lrt <- try(lavaan::lavTestLRT(baseline$model, alternative$model), silent = TRUE)
      if (any(class(lrt)=='try.error')) {
        diffs[paste0(c('delta.chisq.', 'delta.df.', 'delta.pvalue.'), type)] <- NA
      } else {
        tmp <- c(lrt[2, 5:7])
        diffs[paste0(c('delta.chisq.', 'delta.df.', 'delta.pvalue.'), type)] <- tmp
      }
    }
    
    if (class(baseline$model)[1] == 'mplus.model') {
      correction <- ifelse(is.null(baseline$ChiSqM_ScalingCorrection), 1, 
        (baseline$npar * baseline$ChiSqM_ScalingCorrection - alternative$npar*alternative$ChiSqM_ScalingCorrection)/(baseline$npar - alternative$npar))
      diff_chisq <- -2*(baseline$LL - alternative$LL)/correction
      diff_df <- alternative$npar - baseline$npar
      diff_pvalue <- stats::pchisq(diff_chisq, diff_df, lower.tail = FALSE)
      
      diffs[paste0(c('delta.chisq.', 'delta.df.', 'delta.pvalue.'), type)] <- c(diff_chisq, diff_df, diff_pvalue)
    }
  }

  return(diffs)
  
}