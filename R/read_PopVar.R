#' Prepare phenotype and marker genotype data
#' 
#' @description
#' Combines marker genotype data and phenotype data into a single object for downstream
#' analysis and processing.
#' 
#' @param geno An object of class \code{PopVar.geno} from the \code{\link{read_geno}} function.
#' @param pheno A data frame or filename of the phenotypic data. The first column of the data frame or file should
#' be the genotype ID, columns 2 through \code{n.traits + 1} are trait columns, and subsequent columns are fixed effects.
#' @param n.traits The number of traits in the \code{pheno} data frame or file; each trait should be contained in a separate column.
#' @param delim The delimiter character in the pheno data file (if a filename is provided). For example,
#' "," is for csv files, and "\t" is for tab-delimited files.
#' 
#' @details
#' Genotypes in the \code{pheno} file or data frame that are not genotyped are removed from the analysis.
#' 
#' 
#' 
#' @return A variable of class \code{PopVar.data}
#' 
#' @examples
#' # Read in phased cranberry marker data
#' geno_in <- read_geno(geno = cranberry_geno_phased, phased = TRUE, sep = "_hap", min.mac = 5)
#' 
#' # Create the PopVar.data object
#' popvar_data <- read_PopVar(geno = geno_in, pheno = cranberry_pheno, n.traits = 6)
#' 
#' 
#' @export
#' 
read_PopVar <- function(geno, pheno, n.traits, delim = ",") {
  
  # Error checking
  stopifnot(inherits(geno, "PopVar.geno"))
  stopifnot(is.numeric(n.traits))
  n.traits <- as.integer(n.traits)
  stopifnot(is.character(delim))
  
  if (!is.data.frame(pheno)) {
    # Read the file in
    stopifnot(file.exists(pheno))
    
    pheno <- read.table(file = pheno, header = TRUE, sep = delim, as.is = TRUE, check.names = FALSE)
    
  } else {
    pheno <- as.data.frame(pheno)
    
  }
  
  # Match genotypes with phenotypic data
  id.pheno <- unique(pheno[,1])
  id.geno <- colnames(geno@geno.mat)
  id.pheno.geno <- intersect(id.pheno, id.geno)
  
  # Subset the phenotypic data for ids with genotype data
  pheno1 <- pheno[pheno[,1] %in% id.pheno.geno, , drop = FALSE]
  n.id.pheno.geno <- length(id.pheno.geno)
  cat(paste("N =", n.id.pheno.geno, "individuals with phenotypic and genotypic information \n"))
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
  
  # Create the PopVar.data object
  output <- new(Class = "PopVar.data", geno = geno, pheno = pheno2, fixed = fixed, traits = traits)
  return(output)
  
}