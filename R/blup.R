blup1 <- function(phenotype         = phenotype,
                  it                = it,
                  trait_names       = trait_names){
  table_trait <- table(trait_names)
  id <- which(trait_names == names(table_trait)[it])
  phe1 <- t(phenotype[id, , drop = FALSE])
  id_name <- rownames(phe1)
  del_ind_id <- which(apply(phe1, 1, function(x) {
    all(is.na(x))
  }))
  del_env_id <- which(apply(phe1, 2, function(x) {
    all(is.na(x))
  }))
  if (length(del_ind_id) != 0)
    phe1 <- phe1[-del_ind_id, , drop = FALSE]
  if (length(del_env_id) != 0)
    phe1 <- phe1[, -del_env_id, drop = FALSE]
  n <- nrow(phe1)
  ne <- ncol(phe1)
  ind_name <- rownames(phe1)
  env_name <- colnames(phe1)
  if (is.null(env_name))
    env_name <- 1:ne
  phe <- as.numeric(phe1)
  env <- rep(env_name, each = n)
  ind <- rep(ind_name, ne)
  blup <- lme4::lmer(phe ~ (1 | env) + (1 | ind))
  blup_value <- stats::coef(blup)$ind[, 1]
  blup_m <- matrix(blup_value, ncol = 1)
  colnames(blup_m) <- "blup"
  rownames(blup_m) <- ind_name
  blup <- matrix(NA, nrow = length(id_name), ncol = 1)
  colnames(blup) <- "BLUP"
  rownames(blup) <- id_name
  blup[match(rownames(blup_m), rownames(blup)), 1] <- blup_m[, 1]
  return(blup = blup)
}



calc_blup <- function(phenotype         = phenotype,
                      cores             = cores,
                      trait_names       = trait_names) {
  phenotype = phenotype
  cores = cores
  table_trait <- table(trait_names)
  blup_id <- which(table_trait > 1)
  if (length(blup_id) > 0) {
    by_one <- function(it){
      blup1(phenotype         = phenotype,
            it                = it,
            trait_names       = trait_names)
    }
    cores <- setup_cluster(cores = cores, quiet = TRUE)
    blup0 <- cluster_lapply(cores, blup_id, by_one)
    blup0 <-  t(as.matrix(do.call("cbind", blup0)))
    rownames(blup0) <- names(table_trait)[blup_id]
    return(blup0)
  } else {
    return(NULL)
  }
}
