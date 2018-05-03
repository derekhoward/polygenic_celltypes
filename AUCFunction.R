# from https://github.com/sarbal/EGAD/blob/master/EGAD/R/auroc_analytic.R
# by Sara Ballouz
auroc_analytic <- function(scores, labels) {
  
  negatives <- which(labels == 0, arr.ind = TRUE)
  scores[negatives] <- 0
  
  p <- sum(scores, na.rm = TRUE)
  nL <- length(labels)
  np <- sum(labels, na.rm = TRUE)
  nn <- nL - np
  
  auroc <- (p/np - (np + 1)/2)/nn
  
  return(auroc)
} 
