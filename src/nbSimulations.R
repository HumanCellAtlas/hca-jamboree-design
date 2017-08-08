# assume we have the following simulation parameters:

#' @param means vector of means of counts (one for each gene)
#' @param variances vector of variances of counts (one for each gene)
#' @param nSamples scalar that specifies how many samples to generate counts for
#' @param cvRatio positive scalar value indicating a value to multiply the coefficient
#'   of variation
#' @return a matrix of integer counts with genes in rows and samples in columns
#' 
simCountsNB <- function(means, variances, nSamples=100, cvRatio=NULL, ...){
  
  if (is.null(cvRatio)){
    cvRatio <- 1
  }
  
  calcRP <- function (Emean, Evar) 
  {
    RP <- rep(NA, 2)
    if (Evar*cvRatio > Emean){
      # nb parameters
      RP[1] = (Emean)^2/(Evar*cvRatio - Emean)
      RP[2] = Emean/(Evar*cvRatio)
    }else{
      # variance doesn't seem to be greater than the mean; simulate poisson
      RP[1] = Emean
    }
    return(RP)
  }
  
  # put means and variance vectors into a matrix with two columns
  MV <- cbind(means, variances)

  # estimate parameters of negative binomial for each gene
  RP <- t(apply(MV, 1, function(x) calcRP(Emean=x[1], Evar=x[2])))
  
  # given the mean and variance parameters, draw a negative binomial count for each 
  counts <- t(apply(RP, 1, function(x) {
    if (!is.na(x[2])){
      rnbinom(n=nSamples, size=x[1], prob=x[2])
    }else{
      rpois(n=nSamples, lambda=x[1])
    }
  }))
  
  # return a sample x gene matrix of counts
  return(counts)
}

