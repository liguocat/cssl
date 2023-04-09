tidy_data <- function(genotype  = genotype,
                      phenotype = phenotype,
                      dire      = dire,
                      cores     = cores,
                      BLUP      = BLUP)
{
  # cores
  if (is.null(cores) || cores < 1 || cores > parallel::detectCores()) {
    cores <- parallel::detectCores() - 1
    if (cores == 0) cores <- 1
    if (cores > 10) cores <- 10
  }
  cores <- as.integer(cores)
  # A directory for saving results
  if (!is.null(dire)) {
    if (!dir.exists(dire)) {
      suppressWarnings(dir.create(dire, recursive = TRUE))
      if (!dir.exists(dire)) {
        dir1 <- getwd()
      } else {
        dir1 <- dire
      }
    } else{
      dir1 <- dire
    }
  } else {
    dir1 <- getwd()
  }
  if (utils::tail(strsplit(dir1, "")[[1]], 1) != "/")
    dir1 <- paste0(dir1, "/")
  # genotype
  if (is.data.frame(genotype) | is.character(genotype)) {
    if (is.character(genotype)) {
      if (file.exists(genotype)) {
        if (utils::tail(unlist(strsplit(basename(genotype), "[.]")), 1) ==
            "csv") {
          genotype <- data.table::fread(file = genotype, header = TRUE,
                                        sep = ",", check.names = TRUE)
        } else {
          stop("The 'genotype' file must be a comma separated values (csv) file.")
        }

      }
    }
    genotype <- as.data.frame(genotype)
    if (ncol(genotype) == 4 | ncol(genotype) == 5) {
      if (ncol(genotype) == 4) genotype$Type <- "ho"
      names(genotype) <- c("Line", "Chromosome", "Start", "End", "Type")
    } else {
      stop("The 'genotype' file is error.")
    }
    genotype$Type <- tolower(substr(genotype$Type, 1, 2))
    ind_name <- unique(genotype$Line)
    # pmap
    pmap1 <- unique(genotype[, 2:4])
    pmap1 <- pmap1[order(pmap1$Chromosome, pmap1$Start, pmap1$End), ]
    chr <- unique(pmap1$Chromosome)
    # new pmap
    pmap <- NULL
    for (i in 1:length(chr)) {
      # i = 1
      pmap_i <- NULL
      pmap1_i <- pmap1[pmap1$Chromosome == chr[[i]], ]
      pos_i <- sort(unique(c(pmap1_i$Start, pmap1_i$End)))
      for (j in 1:nrow(pmap1_i)) {
        # j = 1
        pmap1_i_j <- pmap1_i[j, ]
        id_ij <- which((pos_i > pmap1_i_j$Start & pos_i <
                          pmap1_i_j$End) == TRUE)
        if (length(id_ij) != 0) {
          pos_pos <- sort(unique(c(pmap1_i_j$Start, pmap1_i_j$End,
                                   pos_i[id_ij])))
          pmap_i_j <- data.frame(Chromosome = c(pmap1_i_j$Chromosome),
                                 Start = pos_pos[-length(pos_pos)],
                                 End = pos_pos[-1])
        }else{
          pmap_i_j <- pmap1_i_j
        }
        pmap_i <- rbind(pmap_i, pmap_i_j)
      }
      pmap <- rbind(pmap, pmap_i)
    }
    pmap <- unique(pmap)
    pmap <- pmap[order(pmap$Chromosome, pmap$Start, pmap$End), ]
    # genotype
    # ind <- data.frame(table(genotype$Line))
    geno <- matrix(-1, nrow = nrow(pmap), ncol = length(ind_name))
    pmap$Marker <- paste0("M", 1:nrow(pmap))
    colnames(geno) <- ind_name
    rownames(geno) <- pmap$Marker
    ID <- 1:nrow(pmap)
    for (i in 1:length(ind_name)) {
      # i = 12
      genotype_i <- genotype[genotype$Line == ind_name[i], ]
      hh <- NULL
      id <- NULL
      for (j in 1:nrow(genotype_i)) {
        # j = 1
        pmap_j <- pmap[pmap$Chromosome == genotype_i$Chromosome[[j]], ]
        ID_j <- ID[pmap$Chromosome == genotype_i$Chromosome[[j]]]
        id_j <- NULL
        for (k in 1:nrow(pmap_j)) {
          # k = 1
          pmap_j_k <- pmap_j[k, ]
          if ((pmap_j_k$Start >= genotype_i$Start[[j]]) &
              (genotype_i$End[[j]] >= pmap_j_k$End)) {
            id_j <- c(id_j, ID_j[k])
          }
        }
        hh <- c(hh, rep(genotype_i$Type[[j]], length(id_j)))
        id <- c(id, id_j)
      }
      if ("ho" %in% hh) hh[hh == "ho"] <- 1
      if ("he" %in% hh) hh[hh == "he"] <- 0
      geno[id, i] <- as.numeric(hh)
    }
  } else if (is.matrix(genotype)) {
    if (is.null(colnames(genotype)) | is.null(rownames(genotype))) {
      stop(paste("The matrix of genotype must have row and column names, the",
                 "row names is marker names and the column names is line",
                 "names."))
    } else {
      geno <- apply(genotype, 2, as.numeric)
      if (!all(unique(as.vector(geno)) %in% c(-1, 0, 1))) {
        stop("The matrix of genotype can only contian -1, 1 and 0.")
      }
    }
    pmap <- data.frame(Marker = colnames(geno))
  } else {
    stop("The argument 'genotype' is error.")
  }
  # phenotype
  if (is.character(phenotype)) {

  }
  if (is.character(phenotype)) {
    if (file.exists(phenotype)) {
      if (utils::tail(unlist(strsplit(basename(phenotype), "[.]")), 1) ==
          "csv") {
        phenotype <- data.table::fread(file = phenotype, header = TRUE,
                                       sep = ",", check.names = TRUE)
        names(phenotype)[1:2] <- c("Trait", "Environment")
      } else {
        stop("The 'phenotype' file must be a comma separated values (csv) file.")
      }
    } else {
      stop("The argument 'phenotype' is error.")
    }
  } else if (is.data.frame(phenotype)) {
    names(phenotype)[1:2] <- c("Trait", "Environment")
  } else {
    stop("The argument 'phenotype' is error.")
  }
  # tidy phenotype
  phenotype <- as.data.frame(phenotype)
  if (nrow(phenotype) == 0) {
    stop("The argument 'phenotype' is error.")
  } else {
    trait_names <- phenotype$Trait
    environment_names <- phenotype$Environment
    phenotype <- as.matrix(phenotype[, -(1:2), drop = FALSE])
    name_col <- colnames(phenotype)
    if (nrow(phenotype) == 1) {
      phenotype <- matrix(as.numeric(phenotype), nrow = 1)
    } else {
      phenotype <- apply(phenotype, 2, function(x){
        suppressWarnings(as.numeric(x))
      })
    }
    colnames(phenotype) <- name_col
    if (BLUP) {
      blup_value <- calc_blup(phenotype         = phenotype,
                              cores             = cores,
                              trait_names       = trait_names)
      if (!is.null(blup_value)) {
        phenotype <- rbind(phenotype, blup_value)
        trait_names <- c(trait_names, rownames(blup_value))
        environment_names <- c(environment_names, rep("BLUP", nrow(blup_value)))
      }
    }
    order1 <- order(trait_names)
    phenotype <- phenotype[order1, , drop = FALSE]
    trait_names <- trait_names[order1]
    environment_names <- environment_names[order1]
  }
  ind_names <- intersect(colnames(geno), colnames(phenotype))
  if (length(ind_names) == 0) {
    stop("The arguments 'genotype' and 'phenotype' are error.")
  }
  geno <- geno[, match(ind_names, colnames(geno)), drop = FALSE]
  phenotype <- phenotype[, match(ind_names, colnames(phenotype)), drop = FALSE]
  return(list(geno              = geno,
              phenotype         = phenotype,
              pmap              = pmap,
              cores             = cores,
              dir1              = dir1,
              trait_names       = trait_names,
              environment_names = environment_names))
}
