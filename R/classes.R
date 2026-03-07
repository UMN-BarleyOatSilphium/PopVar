#' S4 class of genotypic and phenotypic data
#'
#' @slot ploidy The ploidy level.
#' @slot map A data frame of the genetic map information.
#' @slot geno.mat A matrix of marker genotypes coded as allele dosages z = { 0, 1, 2 }
#' @slot allele.freq A numeric vector of allele frequencies of the first allele.
#' @slot haplo.mat A matrix of haplotypes.
#' @slot phased Whether the marker data was phased.
#' @slot missing Whether there is any missing marker data.
#' @slot assume.inbred Whether to assume unphased marker data comes from completely homozygous individuals.
#' @slot dominance Whether to compute dominance marker matrices and effects.
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
                                  allele.freq = "numeric",
                                  haplo.mat = "matrix",
                                  phased = "logical",
                                  missing = "logical",
                                  assume.inbred = "logical",
                                  dominance = "logical",
                                  pheno = "data.frame",
                                  fixed = "data.frame",
                                  traits = "character"))


#' 
#' S4 class of genotypic and phenotypic data with estimated marker effects
#' 
#' @slot marker.effects A list of  marker effect arrays.
#' 
PopVar.me <- setClass(Class = "PopVar.me", 
                      slots = c(marker.effects = "list"),
                      contains = "PopVar.data")

