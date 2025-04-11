#' S4 class of prepared marker genotype data
#'
#' @slot ploidy The ploidy level.
#' @slot map A data frame of the genetic map information.
#' @slot geno.mat A matrix of marker genotypes.
#' @slot haplo.mat A matrix of haplotypes.
#' @slot phased Whether the marker data was phased.
#' @slot G The additive genomic relationship matrix.
#' @slot D The dominance relationship matrix.
#'
PopVar.geno <- setClass(Class = "PopVar.geno", 
                        slots = c(ploidy = "integer", 
                                  map = "data.frame",
                                  geno.mat = "matrix",
                                  haplo.mat = "matrix",
                                  phased = "logical",
                                  G = "Matrix",
                                  D = "Matrix"))


#' S4 class of genotypic and phenotypic data
#'
#' @slot geno A \code{\link{PopVar.geno}} object.
#' @slot pheno A data frame of phenotypic data
#' @slot fixed A data frame of fixed effects
#' @slot traits A list of trait names
#' 
#'
PopVar.data <- setClass(Class = "PopVar.data", 
                        slots = c(geno = "PopVar.geno",
                                  pheno = "data.frame",
                                  fixed = "data.frame",
                                  traits = "character"))


