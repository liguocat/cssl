#' @title Quantitative trait loci (QTLs) mapping for quantitative trait in
#' chromosome segment substitution lines (CSSLs)
#'
#' @description The multi-locus linear model for quantitative traits in CSSLs
#' is described as follow
#' \deqn{y_i = \mu + \sum_{j=1}^m{\omega_{ij}a_j} + \sum_{j=1}^m{\nu_{ij}d_j} +
#'  \varepsilon_i} where the \eqn{y_i} is the phenotypic observation values
#'  of quantitative trait of the *i*th line (\eqn{i = 1, 2, \cdots, n}); \eqn{\mu}
#'  is total average; \eqn{a_j} and \eqn{d_j} are the additive and dominance
#'  effects of the *j*th locus, \eqn{\omega_{ij}} and \eqn{\nu_{ij}} are
#'  their corresponding indicator variables. We assumed that the allele A and a
#'  of each locus are derived from donor and background parents, so
#'  \deqn{\omega_{ij}=\begin{cases}1 & AA \\ 0 & Aa \\ -1 & aa \end{cases}}
#'  and \deqn{\nu_{ij}=\begin{cases}0 & AA \\ 1 & Aa \\ 0 & aa \end{cases}}
#'
#'
#' @param genotype The parameter allows three file input formats.
#' - Character string with path of the genotype file;
#' - Data frame of the genotype;
#' - Matrix of the genotype.
#'
#' See details in vignette.
#' @param phenotype The parameter allows three file input formats.
#' - Character string with path of the phenotype file;
#' - Data frame of the phenotype;
#'
#' See details in vignette.
#' @param dire Character string with path for saving results. Defaults value is
#' \code{NULL} and will save in working directory, using [getwd] to get the
#' working directory.
#' @param cores The number of CPU cores to use, for parallel calculations.
#' Defaults is 1 and no parallel.
#' @param LOD A numerical value indicates the threshold of LOD. Defaults is
#' \code{2.5}.
#' @param BLUP Defaults is \code{FALSE}. If it is \code{TRUE}, will calculate
#' the BLUP (best linear unbiased prediction) values for trait measured in
#' multi-environment.
#' @param qq Defaults is \code{FALSE}. If it is \code{TRUE}, will transform
#' quantitative trait phenotypic value using function [stats::qqnorm]. It is
#' useful for expression quantitative trait locus (eQTL) mapping.
#'
#' @return A list include two data frame.
#'
#' \item{map}{The map have four columns, including:}
#' - Chromosome, the names of chromosome
#' - Start, the position of start
#' - End, the position of end
#' - Marker, the names of marker
#'
#' \item{result}{The result have eight columns, including:}
#' - Trait, the trait names
#' - Environment, the environment names
#' - Marker, the marker names
#' - Effect, the values of effect
#' - Type, "additive" or "dominant"
#' - LOD, the LOD value
#' - Same_marker, the code of Same_marker is same as the Marker
#' - r^2(%), the ratio of QTL variance and total variance
#'
#' @export
#'
#' @examples
#' # library(cssl)
#' \dontrun{
#' #' data(SimData)
#' result <- cssl(genotype  = SimData$genotype,
#'                phenotype = SimData$phenotype,
#'                LOD       = 5)
#' }

cssl <- function(genotype  = genotype,
                 phenotype = phenotype,
                 dire      = NULL,
                 cores     = 1,
                 LOD       = 2.5,
                 BLUP      = FALSE,
                 qq        = FALSE)
{
  data_all <- tidy_data(genotype  = genotype,
                        phenotype = phenotype,
                        dire      = dire,
                        cores     = cores,
                        BLUP      = BLUP)
  geno <- data_all$geno
  pmap <- data_all$pmap
  phenotype <- data_all$phenotype
  trait_names <- data_all$trait_names
  environment_names <- data_all$environment_names
  cores <- data_all$cores
  dir1 <- data_all$dir1
  geno_a <- geno
  geno_d <- 1 - abs(geno)
  result <- mapping(phenotype         = phenotype,
                    geno_a            = geno_a,
                    geno_d            = geno_d,
                    trait_names       = trait_names,
                    environment_names = environment_names,
                    qq                = qq,
                    LOD               = LOD,
                    cores             = cores)
  result <- list(map = pmap, result = result)
  it <- 1
  while (file.exists(paste0(dir1, "result_", it, ".xlsx"))) {
    it <- it + 1
  }
  openxlsx::write.xlsx(result, paste0(dir1, "result_", it, ".xlsx"),
                       overwrite = TRUE, colName = TRUE)
  base::cat("The result saved in ", paste0(dir1, "result_", it, ".xlsx"), "\n")
  return(result)
}
