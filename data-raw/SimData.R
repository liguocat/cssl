## code to prepare `SimData` dataset goes here

genotype <- read.table("D:/CSSL/genotype.csv", sep = ",", header = TRUE)
genotype$Chromosome <- match(genotype$Chromosome, unique(genotype$Chromosome))
genotype <- genotype[genotype$Chromosome %in%
                       sample(x = unique(genotype$Chromosome), size = 12), ]
genotype$Chromosome <- match(genotype$Chromosome, unique(genotype$Chromosome))
genotype$Type[genotype$Type == "homo"] <- "homozygote"
genotype$Type[genotype$Type == "hetero"] <- "heterozygote"
genotype$Line <- as.numeric(substr(genotype$Line, 2, 10))
genotype <- genotype[genotype$Line %in% sample(x = unique(genotype$Line),
                                               size = 100), ]
genotype$Line <- match(genotype$Line, unique(genotype$Line))
genotype <- genotype[order(genotype$Line, genotype$Chromosome,
                           genotype$Start), ]
genotype$Line <- paste0("cssl", genotype$Line)
if (!dir.exists("inst/extdata/")) {
  dir.create("inst/extdata/", recursive = TRUE)
}
write.table(genotype, "inst/extdata/CSSL100.csv", col.names = TRUE,
            row.names = FALSE, sep = ",")

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

geno_a <- geno
geno_a <- unique(geno_a)
geno_d <- 1 - abs(geno)
geno_d <- unique(geno_d)
del_id <- apply(geno_d, 1, function(x){all(x == 0)})
if (length(del_id) != 0) geno_d <- geno_d[-del_id, ]

eff_a <- NULL
eff_d <- NULL
id_a <- c(33, 67)
id_d <- c(5, 12)
for (i in 1:2) {
  fn <- function(x){abs(var(geno_a[id_a[i], ] * x) - 1)}
  eff_a <- c(eff_a, optim(par = 1, fn = fn, method = "L-BFGS-B",
                          lower = -15, upper = 15)$par * (-1)^i)
}
for (i in 1:2) {
  fn <- function(x){abs(var(geno_d[id_d[i], ] * x) - 1)}
  eff_d <- c(eff_d, optim(par = 1, fn = fn, method = "L-BFGS-B",
                          lower = -15, upper = 15)$par * (-1)^i)
}

y <- t(geno_a[c(33, 67), ]) %*% matrix(eff_a) +
  t(geno_d[c(5, 12), ]) %*% matrix(eff_d)

phenotype <- NULL
for (i in 1:7) {
  phenotype <- rbind(phenotype, matrix(y + rnorm(100), nrow = 1))
}
phenotype <- data.frame(Trait = paste0("Trair", c(rep(1, 3), rep(2, 4))),
                        Environment = paste0("Env", c(1:3, 1:4)),
                        phenotype)
names(phenotype) <- c("Trait", "Environment", paste0("cssl", 1:100))
write.table(phenotype, "inst/extdata/phenotype.csv", col.names = TRUE,
            row.names = FALSE, sep = ",")

SimData <- list(genotype  = genotype,
                phenotype = phenotype,
                geno      = geno)

usethis::use_data(SimData, overwrite = TRUE)
