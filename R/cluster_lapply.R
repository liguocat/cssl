# Reference from R/qtl2 package
# run code by cluster (generalizes lapply, parLapply, and mclapply)
# (to deal with different methods on different architectures)
# if cores==1, just use lapply
setup_cluster <- function(cores, quiet = TRUE) {
  if (cores > 1 && Sys.info()[1] == "Windows") {
    # windows doesn't support mclapply
    cores <- parallel::makeCluster(cores)
    # the following calls on.exit() in the function that called this one
    # see http://stackoverflow.com/a/20998531
    do.call("on.exit",
            list(quote(parallel::stopCluster(cores))),
            envir = parent.frame())
  }
  return(cores)
}




cluster_lapply <- function(cores, ...){
  if (!is.numeric(cores)) {
    # cluster object; use mclapply
    return(parallel::parLapply(cores, ...))
  } else {
    if (cores == 1) return(lapply(...))
    return(parallel::mclapply(..., mc.cores = cores))
  }
}
