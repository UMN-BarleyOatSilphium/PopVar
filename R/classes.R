#' S4 class of genotypic and phenotypic data
#'
#' @slot ploidy The ploidy level.
#' @slot map A data frame of the genetic map information.
#' @slot geno.mat A matrix of marker genotypes.
#' @slot haplo.mat A matrix of haplotypes.
#' @slot phased Whether the marker data was phased.
#' @slot missing Whether there is any missing marker data.
#' @slot assume.inbred Whether to assume unphased marker data comes from completely homozygous individuals.
#' @slot geno A \code{\link{PopVar.geno}} object.
#' @slot pheno A data frame of phenotypic data
#' @slot fixed A data frame of fixed effects
#' @slot traits A list of trait names
#' 
#' @rdname PopVar.data-class
#'
PopVar.data <- setClass(Class = "PopVar.data", 
                        slots = c(ploidy = "integer", 
                                  map = "data.frame",
                                  geno.mat = "matrix",
                                  haplo.mat = "matrix",
                                  phased = "logical",
                                  missing = "logical",
                                  assume.inbred = "logical",
                                  pheno = "data.frame",
                                  fixed = "data.frame",
                                  traits = "character"))


#' 
#' S4 class of genotypic and phenotypic data with estimated marker effects
#' 
#' @slot marker.effects A matrix of marker effects
#' 
PopVar.me <- setClass(Class = "PopVar.me", 
                      slots = c(marker.effects = "array"),
                      contains = "PopVar.data")

