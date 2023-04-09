LogLikelihood <- function(y, expec, sigma2){
  n <- length(y)
  LogL <- -n * log(sqrt(2 * pi)) -
    0.5 * n * log(sigma2) - 0.5 * crossprod(y - expec) / sigma2
  return(as.numeric(LogL))
}
