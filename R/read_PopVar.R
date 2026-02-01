#' Read in phenotypic and marker genotype data
#' 
#' @description
#' Reads and combines marker genotype data and phenotype data into a single object for downstream
#' analysis and processing.
#' 
#' @param geno A data frame or path to a file containing the marker allele dosages. See Details for information
#' about formatting.
#' @param pheno A data frame or path to a file containing the phenotypic data. See Details for information
#' about formatting.
#' @param n.traits The number of traits in the \code{pheno} data frame or file; each trait should be contained in a separate column.
#' @param ploidy The ploidy level e.g. 2, 4, 6, ..., etc.
#' @param phased Logical. If TRUE, the input marker genotype data is assumed phased, with alleles coded as positive integers. 
#' The name of each individual/clone column is the id and haplotype, concatenated by \code{sep}.
#' @param sep The haplotype concatenator; see \code{phased}.
#' @param assume.inbred Logical. If TRUE, individuals/clones are assumed to be completely homozygous. Setting \code{phased = FALSE}
#' \code{assume.inbred = FALSE} is useful for testing cross-validation, but you will not be able to predict genetic variances. 
#' This argument is ignored if \code{phased = TRUE}.
#' @param min.mac The minimum minor allele count to retain a marker.
#' @param max.ind.miss The maximum rate of missingness to retain an individual.
#' @param max.mar.miss The maximum rate of missingness to retain an marker.
#' @param max.r2 The maximum squared correlation between any pair of markers. The \code{r^2} is used directly in 
#' filtering. For each pair of markers exceeding this threshold, the one with the greater minor allele frequency is 
#' retained. Set \code{max.r2 = 1} for no filtering.
#' @param dominance Logical. Should a dominance relationship matrix be constructed?
#' @param delim The delimiter character for the \code{geno} or \code{pheno} data files (if filenames are provided). For example,
#' "," is for csv files, and "\t" is for tab-delimited files.
#' 
#' 
#' @details
#' 
#' The first columns of the file or \code{geno} data frame or file should be marker, chrom, cM (and optionally bp, if provided). Subsequent columns 
#' contain the allele dosage for individuals/clones/haplotypes, coded 0,1,2, ..., \code{ploidy}. While fractional 
#' values are allowed for genomewide prediction, these values are rounded when predicting the genetic variance within crosses.
#' 
#' If unphased marker data is passed, the function will assume that all individuals are homozygous at all sites. Marker
#' genotype doses will be rounded to meet that constraint.
#' 
#' Missing allele dosage information is not permitted for downstream analysis. If there are missing allele dosages,
#' the user will be prompted to run genotype imputation.
#' 
#' The first column of the \code{pheno} data frame or file should be genotype ID, columns 2 through \code{n.traits + 1} 
#' are trait columns, and subsequent columns are fixed effects.
#' 
#' Individuals/clones in the \code{pheno} file or data frame that are not genotyped are removed from the analysis.
#' 
#' Only diploid allele dosages are allowed at this time.
#' 
#' When \code{dominance=FALSE}, non-additive effects are not captured.
#' If \code{dominance=TRUE}, a (digenic) dominance covariance matrix is used instead.
#' 
#' 
#' @return A variable of class \code{\link[PopVar.data]{PopVar.data}}
#' 
#' @examples
#' # Read in phased cranberry marker data and phenotypes
#' geno <- cranberry_geno_phased
#' pheno <- cranberry_pheno
#' dat <- read_PopVar(geno = geno, pheno = pheno, n.traits = 6, phased = TRUE, sep = "_")
#' 
#' # Read in unphased barley marker data and phenotypes; this will assume that all
#' # individuals are completely homozygous
#' geno <- barley_geno
#' pheno <- barley_pheno
#' dat <- read_PopVar(geno = geno, pheno = pheno, n.traits = 4, phased = FALSE, assume.inbred = TRUE)
#' 
#' # Read in unphased cranberry marker data and phenotypes
#' geno <- cranberry_geno
#' pheno <- cranberry_pheno
#' dat <- read_PopVar(geno = geno, pheno = pheno, n.traits = 6, phased = FALSE, assume.inbred = FALSE)
#' 
#' @export
#' 
read_PopVar <- function(geno, pheno, n.traits, ploidy = 2L, phased = FALSE, sep = ".", assume.inbred = FALSE, 
                        max.ind.miss = 0.5, max.mar.miss = 0.5, min.mac = 5, max.r2 = 1, 
                        dominance = FALSE, delim = ",") {
  
  # Error checking
  ploidy <- as.integer(ploidy)
  stopifnot(is.logical(phased))
  stopifnot(is.logical(dominance))
  stopifnot(is.character(sep))
  stopifnot(is.numeric(n.traits))
  n.traits <- as.integer(n.traits)
  stopifnot(is.character(delim))
  
  if (ploidy > 2) stop("PopVar only supports diploids currently.")
  
  # Check geno
  if (!is.data.frame(geno)) {
    # Read the file in
    stopifnot(file.exists(geno))
    
    geno <- read.table(file = geno, header = TRUE, sep = delim, as.is = TRUE, check.names = FALSE)
    
  } else {
    geno <- as.data.frame(geno)
    
  }
  
  
  # Check pheno
  if (!is.data.frame(pheno)) {
    # Read the file in
    stopifnot(file.exists(pheno))
    
    pheno <- read.table(file = pheno, header = TRUE, sep = delim, as.is = TRUE, check.names = FALSE)
    
  } else {
    pheno <- as.data.frame(pheno)
    
  }
  

  # Extract the map information
  cols_required <- c("marker", "chrom", "cM")
  if (any(!cols_required %in% colnames(geno))) {
    stop(paste0("The first three columns of 'geno' must be '", paste0(cols_required, collapse = "', '"), "'."))
  }
  if (colnames(geno)[4] == "bp") {
    map <- geno[,1:4]
    geno <- as.matrix(geno[,-(1:4)])
  } else {
    map <- geno[,1:3]
    geno <- as.matrix(geno[,-(1:3)])
  }
  row.names(geno) <- map$marker
  
  # Detect any missing data
  if (any(is.na(geno))) {
    missing <- TRUE
  } else {
    missing <- FALSE
  }
  
  # Keep a copy of the original genos and map
  geno_mat_orig <- geno
  map_orig <- map
  
  # Make a haplotype matrix if not phased; make a genotype matrix if phase
  if (!phased) {
    geno_mat <- geno
    
    # If assume inbred is TRUE, eliminate hets in the genotype matrix and create a haplotype matrix
    if (assume.inbred) {
      geno_mat <- ifelse(geno_mat >= (ploidy / 2), ploidy, 0)
    }
    
    # Expand to a haplotype matrix
    haplo_mat <- matrix(as.numeric(NA), nrow = nrow(geno_mat), ncol = ncol(geno_mat) * 2)
    
    hap_prefix <- "hap"
    hap_prefixes <- paste0("_", hap_prefix, seq_len(ploidy))
    
    id <- colnames(geno_mat)
    colnames(haplo_mat) <- paste0(rep(id, each = 2), rep(hap_prefixes, length(id)))
    # Iterate over id
    for (idx in id) {
      haplo_idx <- cbind(geno_mat[,idx], geno_mat[,idx])
      haplo_idx <- haplo_idx / 2
      colnames(haplo_idx) <- paste0(idx, hap_prefixes)
      haplo_mat[,colnames(haplo_idx)] <- haplo_idx
    }
    
  } else {
    
    # Sum haplotypes to get genotypes
    haplo_mat <- geno
    # Unique column names
    id_names <- colnames(haplo_mat)
    id_names_split <- strsplit(x = id_names, split = sep)
    hap_prefixes <- paste0(sep, unique(sapply(X = id_names_split, FUN = "[[", 2)))
    
    id_names_split <- sapply(X = id_names_split, FUN = "[[", 1)
    id <- unique(id_names_split)
    
    geno_mat <- matrix(data = as.numeric(NA), nrow = nrow(haplo_mat), ncol = length(id), 
                       dimnames = list(row.names(haplo_mat), id))
    
    # Iterate over unique ids
    for (idx in id) {
      # Find the position of the matching haplotypes; sum the haplotype dosages
      geno_mat[,idx] <- rowSums(haplo_mat[,which(id_names_split == idx)])
    }
    
    # Round the haplotype dosages
    haplo_mat <- round(haplo_mat)
    
  }
  
  row.names(haplo_mat) <- row.names(geno_mat) <- row.names(geno)
  
  
  
  # First message
  if (phased) {
    cat("Using phased marker genotypes.\n\n")
  } else {
    if (assume.inbred) {
      cat("Using inferred phased marker genotypes assuming completely inbred genotypes.\n\n")
    } else {
      cat("Using unphased marker genotypes.\n\n")
    }
  }
  
  # Filter
  map_in <- map[map$marker %in% row.names(geno_mat), 1:3]
  colnames(map_in) <- c("marker", "chrom", "position")
  geno_in <- cbind(map_in, geno_mat)
  cat(sub("X", nrow(geno_mat), "Number of markers in the original dataset = X\n"))
  cat(sub("X", ncol(geno_mat), "Number of individuals in the original dataset = X\n"))
  filter_out <- filter_geno(geno = geno_in, max.ind.miss = max.ind.miss, max.mar.miss = max.mar.miss, min.mac = min.mac,
                            max.r2 = max.r2)
  
  idx <- which(row.names(geno_mat) %in% filter_out$marker)
  idy <- which(colnames(geno_mat) %in% names(filter_out))
  
  cat(sub("X", length(idx), "\nNumber of markers in the filtered dataset = X\n"))
  cat(sub("X", length(idy), "Number of individuals in the filtered dataset = X\n"))
  
  # Subset the matrices
  geno_mat <- geno_mat[idx, idy]
  map <- map[idx, ]
  idyy <- apply(X = expand.grid(x = hap_prefixes, y = colnames(geno_mat))[,c(2,1)], MARGIN = 1, FUN = paste0, collapse = "")
  haplo_mat <- haplo_mat[idx, idyy]
  
  
  # Remove genotypes with completely missing phenotypic data
  pheno_orig <- pheno
  idx <- rowMeans(is.na(pheno[,2:(n.traits + 1)])) < 1
  pheno <- pheno[idx, ]
  
  # Match genotypes with phenotypic data
  id.pheno <- unique(pheno[,1])
  id.geno <- colnames(geno_mat)
  id.pheno.geno <- intersect(id.pheno, id.geno)
  
  # Subset the phenotypic data for ids with genotype data
  pheno1 <- pheno[pheno[,1] %in% id.pheno.geno, , drop = FALSE]
  n.id.pheno.geno <- length(id.pheno.geno)
  cat(paste("\nN =", n.id.pheno.geno, "individuals with phenotypic and genotypic information \n"))
  id.geno.only <- setdiff(id.geno, id.pheno)
  n.geno.only <- length(id.geno.only)
  if (n.geno.only > 0) cat(paste("N =", n.geno.only, "individuals with only genotypic information \n"))
  
  # How many individuals with only phenotypic information were removed?
  id.pheno.only <- setdiff(id.pheno, id.geno)
  n.pheno.only <- length(id.pheno.only)
  if (n.pheno.only > 0) cat(paste("N =", n.pheno.only, "individuals with only phenotypic information; these individuals are excluded from the dataset.\n"))
  
  # How many fixed effects?
  n.fixed <- ncol(pheno1) - n.traits - 1
  
  if (n.fixed > 0) {
    fixed <- data.frame(pheno1[,(n.traits+2):ncol(pheno1)], stringsAsFactors = TRUE)
    fixed.names <- colnames(pheno1)[(n.traits+2):ncol(pheno1)]
    colnames(fixed) <- fixed.names
    pheno2 <- data.frame(pheno1[,1:(1+n.traits)], stringsAsFactors = FALSE)
    cat(paste("\nThe following fixed effects were identified in the data:\n",paste(fixed.names,collapse="\n"),"\n",sep=""))
  } else {
    fixed <- data.frame(NULL)
    pheno2 <- pheno1
    cat("\nNo fixed effects were identified in the data.\n")
  }
  traits <- colnames(pheno2)[-1]
  cat(paste("\nThe following traits were identified in the data:\n",paste(traits,collapse="\n"),"\n",sep=""))
  
  # Eliminate the haplotype matrix if assume.inbred is FALSE and phased is FALSE
  if (!phased & !assume.inbred) {
    haplo_mat <- matrix(NA, nrow = 0, ncol = 0)
  }
  
  # Create the PopVar.data object
  output <- new(Class = "PopVar.data", ploidy = as.integer(ploidy), map = map, geno.mat = geno_mat, haplo.mat = haplo_mat, 
                phased = phased, missing = missing, assume.inbred = assume.inbred, pheno = pheno2, fixed = fixed, traits = traits)
  
  return(output)
  
}
  
  
  
  
  