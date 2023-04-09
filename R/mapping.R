mapping1 <- function(phenotype          = phenotype,
                     geno_a             = geno_a,
                     geno_d             = geno_d,
                     trait_names        = trait_names,
                     environment_names  = environment_names,
                     qq                 = qq,
                     LOD                = LOD,
                     it                 = it)
{
  y <- phenotype[it, ]
  tname <- trait_names[it]
  ename <- environment_names[it]
  scan1(y      = y,
        geno_a = geno_a,
        geno_d = geno_d,
        tname  = tname,
        ename  = ename,
        qq     = qq,
        LOD    = LOD)
}





mapping <- function(phenotype          = phenotype,
                    geno_a             = geno_a,
                    geno_d             = geno_d,
                    trait_names        = trait_names,
                    environment_names  = environment_names,
                    qq                 = qq,
                    LOD                = LOD,
                    cores              = cores)
{
  phenotype          = phenotype
  geno_a             = geno_a
  geno_d             = geno_d
  trait_names        = trait_names
  environment_names  = environment_names
  qq                 = qq
  LOD                = LOD
  cores              = cores
  by_one <- function(it){
    mapping1(phenotype          = phenotype,
             geno_a             = geno_a,
             geno_d             = geno_d,
             trait_names        = trait_names,
             environment_names  = environment_names,
             qq                 = qq,
             LOD                = LOD,
             it                 = it)
  }
  cores <- setup_cluster(cores = cores, quiet = TRUE)
  result <- cluster_lapply(cores, 1:nrow(phenotype), by_one)
  result <-  do.call("rbind", result)
  return(result)
}
