#' Read in marker genotype data
#'
#' @param filename Name of CSV file with marker allele dosage
#' @param geno A data frame of marker allele dosages. If passed, \code{filename} is ignored.
#' @param ploidy The ploidy level e.g. 2, 4, 6, ..., etc.
#' @param phased If TRUE, the input genotype matrix is phased, with alleles coded as positive integers. The name of each
#' individual/clone column is the id and haplotype, concatenated by \code{sep}.
#' @param sep The haplotype concatenator; see \code{phased}.
#' @param min.mac The minimum minor allele count to retain a marker.
#' @param max.r The maximum correlation between any pair of markers. The \code{r^2} is used directly in filtering. For each pair of markers exceeding this threshold, the one with the greater minor allele frequency is retained. Set \code{max.r = 1} for no filtering.
#' @param dominance Logical. Should a dominance relationship matrix be constructed?
#' @param inbred Logical. If unphased marker data is passed, should the function assume that all genotypes are inbred?
#' 
#' @details
#' The first columns of the file or geno object are marker, chrom, cM (and optionally bp, if provided).Subsequent columns 
#' contain the allele dosage for individuals/clones, coded 0,1,2. While fractional values are allowed for genomewide prediction, 
#' these values are rounded when predicting the genetic variance within crosses.
#' 
#' Only diploid allele dosages are allowed.
#' 
#' When \code{dominance=FALSE}, non-additive effects are not captured.
#' If \code{dominance=TRUE}, a (digenic) dominance covariance matrix is used instead.
#' 
#' @references Endelman, J.B. 2023. Fully efficient, two-stage analysis of multi-environment trials with directional dominance and multi-trait genomic selection. Theor. Appl. Genet. 136(4): 65.
#' 
#'
#' @return 
#' Variable of class \code{\link{PopVar.geno}}.
#' 
#' @examples
#' # Read in barley marker genotypes assuming no phasing and non-inbred lines.
#' geno_in <- read_geno(geno = barley_geno, min.mac = 5)
#' 
#' # Read in barley marker genotypes assuming completely inbred lines
#' geno_in <- read_geno(geno = barley_geno, min.mac = 5, inbred = TRUE)
#' 
#' # Read in cranberry marker genotype assuming no phasing
#' geno_in <- read_geno(geno = cranberry_geno, min.mac = 5)
#' 
#' # Read in cranberry phased haplotypes from Beagle
#' geno_in <- read_geno(geno = cranberry_geno_phased, phased = TRUE, sep = "_hap", min.mac = 5)
#'
#' @importFrom utils read.csv capture.output
#' @import Matrix
#' @import StageWise
#' 
#' @export
#'
read_geno <- function(filename, geno, ploidy = 2L, phased = FALSE, sep = ".", min.mac = 5, max.r = 1, dominance = FALSE, inbred = FALSE) {
  
  # Error checking
  ploidy <- as.integer(ploidy)
  stopifnot(is.logical(phased))
  stopifnot(is.logical(dominance))
  stopifnot(is.logical(inbred))
  stopifnot(min.mac >= 0)
  stopifnot(max.r >= 0 & max.r <= 1)
  stopifnot(is.character(sep))
  
  
  if (!missing(geno)) {
    data <- geno
  } else {
    data <- read.csv(file = filename, check.names = FALSE)
  }
  
  # Extract the map information
  cols_required <- c("marker", "chrom", "cM")
  if (any(!cols_required %in% colnames(data))) {
    stop("The first three columns of 'filename' or 'geno' must be 'marker', 'chrom', and 'bp'.")
  }
  if (colnames(data)[4] == "bp") {
    map <- data[,1:4]
    geno <- as.matrix(data[,-(1:4)])
  } else {
    map <- data[,1:3]
    geno <- as.matrix(data[,-(1:3)])
  }
  
  m <- nrow(geno)
  
  if (phased) {
    # No missing data is allowed
    if (any(is.na(geno))) stop("Missing data is not allowed with phased marker genotypes.")
    
    # Sum haplotypes to get genotypes
    haplo_mat <- geno
    # Unique column names
    id_names <- colnames(haplo_mat)
    id_names_split <- strsplit(x = id_names, split = sep)
    id_names_split <- sapply(X = id_names_split, FUN = "[[", 1)
    id <- unique(id_names_split)
    
    geno_mat <- matrix(data = as.numeric(NA), nrow = nrow(haplo_mat), ncol = length(id), dimnames = list(row.names(haplo_mat), id))
    
    # Iterate over unique ids
    for (idx in id) {
      # Find the position of the matching haplotypes; sum the haplotype dosages
      geno_mat[,idx] <- rowSums(haplo_mat[,which(id_names_split == idx)])
    }
    
    row.names(haplo_mat) <- row.names(geno_mat) <- data[,1]
    
    # Round the haplotype dosages
    haplo_mat <- round(haplo_mat)
    
  } else if (inbred) {
    # This functionality only works with diploids
    if (ploidy > 2) stop("Assuming inbred lines works only when ploidy = 2.")
    
    # Eliminate hets in the genotype matrix and create a haplotype matrix
    geno_mat <- geno
    # Impute with the population mean
    if (any(is.na(geno_mat))) {
      geno_mat1 <- apply(X = geno_mat, MARGIN = 1, FUN = function(snp) {
        snp_mean <- mean(snp, na.rm = TRUE)
        snp[is.na(snp)] <- snp_mean
        snp
      })
      geno_mat <- t(geno_mat1)
    }
    
    geno_mat <- ifelse(geno_mat >= 1, 2, 0)
    # Expand to a haplotype matrix
    haplo_mat <- matrix(as.numeric(NA), nrow = nrow(geno_mat), ncol = ncol(geno_mat) * 2)
    
    id <- colnames(geno_mat)
    colnames(haplo_mat) <- paste0(rep(id, each = 2), rep(c("_hapA", "_hapB"), length(id)))
    # Iterate over id
    for (idx in id) {
      haplo_idx <- cbind(geno_mat[,idx], geno_mat[,idx])
      haplo_idx <- haplo_idx / 2
      colnames(haplo_idx) <- paste0(idx, c("_hapA", "_hapB"))
      haplo_mat[,colnames(haplo_idx)] <- haplo_idx
    }
    
    row.names(haplo_mat) <- row.names(geno_mat) <- data[,1]
    
  } else {
    haplo_mat <- matrix(data = as.numeric(0), nrow = 0, ncol = 0)
    geno_mat <- geno
    # Impute with the population mean
    if (any(is.na(geno_mat))) {
      geno_mat1 <- apply(X = geno_mat, MARGIN = 1, FUN = function(snp) {
        snp_mean <- mean(snp, na.rm = TRUE)
        snp[is.na(snp)] <- snp_mean
        snp
      })
      geno_mat <- t(geno_mat1)
      geno_mat <- round(geno_mat)
      
    }
    row.names(geno_mat) <- data[,1]
    
  }
  
  id <- colnames(geno_mat)
  p <- rowMeans(x = geno_mat, na.rm = TRUE) / ploidy
  n.minor <- vector("numeric", length = length(p))
  n.minor[is.na(p)] <- as.numeric(NA)
  allele.count <- apply(geno_mat,1,function(z) table(factor(round(z),levels=0:ploidy)) )
  n.minor[which(p > 0.5)] <- colSums(allele.count[, which(p > 0.5)]) - allele.count[ploidy + 1, which(p > 0.5)]
  n.minor[which(p <= 0.5)] <- allele.count[1, which(p <= 0.5)]
  
  ix <- which(n.minor >= min.mac)
  m <- length(ix)
  if (m < 1) stop("No markers passed the 'min.mac' threshold.")
  
  if (phased) {
    cat("Using phased marker genotypes.\n")
  } else if (inbred) {
    cat("Using inferred phased marker genotypes assuming completely inbred genotypes.\n")
  } else {
    cat("Using unphased marker genotypes.\n")
  }
  
  cat(sub("X",min.mac,"Minor allele threshold = X genotypes\n"))
  cat(sub("X",m,"Number of markers passing the 'min.mac' threshold = X\n"))
  geno_mat <- geno_mat[ix,]
  if (nrow(haplo_mat) > 0) {
    haplo_mat <- haplo_mat[ix, ]
  }
  
  p <- p[ix]
  map <- map[ix,]
  maf <- pmin(p, 1-p)
  
  # Calculate the pairwise correlation between markers
  if (max.r < 1) {

    geno_mat1 <- t(geno_mat)
    max.r2 <- max.r^2
    corr <- cor(x = geno_mat1)^2
    corr[lower.tri(corr, diag = TRUE)] <- NA
    
    
    # While loop
    while (any(corr > max.r2, na.rm = TRUE)) {
      # Find the pair with the highest correlation
      pairs <- which(corr == max(corr, na.rm = TRUE), arr.ind = TRUE)
      pair_test <- pairs[1,, drop = FALSE]
      # Find the marker with the highest MAF
      mar_test <- c(row.names(corr)[pair_test[,1]], colnames(corr)[pair_test[,2]])
      remove <- names(which.min(maf[mar_test]))
      
      # Recalculate the correlation
      geno_mat1 <- geno_mat1[, setdiff(colnames(geno_mat1), remove)]
      corr <- cor(x = geno_mat1)^2
      corr[lower.tri(corr, diag = TRUE)] <- NA
      
    }
    
    # Markers to keep
    markers_keep <- colnames(geno_mat1)
    
    ix <- match(markers_keep, row.names(geno_mat))
    p <- p[ix]
    map <- map[ix,]
    geno_mat <- geno_mat[ix,]
    
    m <- length(ix)
    if (m < 1) stop("No markers passed the 'min.mac' threshold.")
    
    cat(sub("X", max.r2, "Maximum marker r2 = X\n"))
    cat(sub("X", m, "Number of markers passing the 'max.r' threshold = X\n"))
    
  }

  
  n <- ncol(geno_mat)
  cat(sub("X",n,"Number of genotypes = X\n"))

  
  coeff <- Matrix::Matrix(scale(t(geno_mat),scale=F), dimnames=list(id,rownames(geno_mat)))
  coeff[is.na(coeff)] <- 0
  
  scale <- ploidy*sum(p*(1-p))
  G <- Matrix::tcrossprod(coeff)/scale
  w <- 1e-5
  
  H <- (1-w)*G + w*mean(diag(G))*Diagonal(n=nrow(G))

  if (dominance) {
    Pmat <- kronecker(matrix(p,nrow=1,ncol=m),matrix(1,ncol=1,nrow=n))
    coeff.D <- Matrix(-2*choose(ploidy,2)*Pmat^2 + 2*(ploidy-1)*Pmat*t(geno_mat) - t(geno_mat)*(t(geno_mat)-1))
    coeff.D[is.na(coeff.D)] <- 0
    scale.D <- 4*choose(ploidy,2)*sum(p^2*(1-p)^2)
    D <- tcrossprod(coeff.D)/scale.D
    D <- (1-w)*D + (w)*mean(diag(D))*Diagonal(n=nrow(D))
    
  } else {
    D <- Diagonal(0)

  }
  
  output <- new(Class = "PopVar.geno", ploidy = as.integer(ploidy), map = map, geno.mat = geno_mat, haplo.mat = haplo_mat, 
                phased = phased, G = G, D = D)
  
  return(output)
}
