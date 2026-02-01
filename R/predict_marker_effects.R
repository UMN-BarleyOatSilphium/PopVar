#' Predict marker effects
#' 
#' 
#' @param data A \code{PopVar.data} object from \code{read_PopVar}.
#' @param models The list of models to run. See Details for supported models.
#' @param params Parameters to pass to \code{\link{BGLR}}, if \code{models} includes Bayesian models.
#' 
#' 
#' @details
#' The function currently supports the models:  
#' \describe{
#'    \item{rrBLUP}{Ridge regression best linear unbiased prediction.}
#'    \item{BayesA}{BayesA.}
#'    \item{BayesB}{BayesB.}
#'    \item{BayesC}{BayesC-pi.}
#'    \item{BL}{Bayesian Lasso.}
#'    \item{BRR}{Bayesian ridge regression.}
#' }
#' 
#' @return A variable of class \code{\link[PopVar.me]{PopVar.me}}
#' 
#' @examples
#' # Read in unphased cranberry marker data and phenotypes
#' geno <- cranberry_geno
#' pheno <- cranberry_pheno
#' dat <- read_PopVar(geno = geno, pheno = pheno, n.traits = 6, phased = FALSE, 
#'                    assume.inbred = FALSE)
#' 
#' # Compute marker effects using the rrBLUP and RKHS models
#' data_me <- predict_marker_effects(data = popvar_data, models = c("rrBLUP", "BRR"))
#' 
#' @importFrom rrBLUP mixed.solve
#' @importFrom BGLR BGLR
#' 
#' @export
#' 
predict_marker_effects <- function(data, models = c("rrBLUP", "BayesA", "BayesB", "BayesC", "BL", "BRR"),
                                   params = list(nIter = 500, burnIn = 50)) {
  
  # Error handling
  stopifnot(inherits(data, "PopVar.data"))
  models <- match.arg(models, several.ok = TRUE)
  stopifnot(is.list(params))
  required_param_names <- c("nIter", "burnIn")
  if (!all(required_param_names %in% names(params))) stop("The following must be in params: ", paste0(required_param_names, collapse = ", "))
  
  # Get the genotype data
  geno <- data@geno.mat
  # Get phenotypes
  pheno <- data@pheno
  pheno$geno_id <- factor(pheno$geno_id, levels = colnames(geno))
  # Get traits
  traits <- data@traits
  
  # Create an array of m x p x r dimensions, where m = number of markers, p is number of traits, and r is number of models
  marker_effects <- array(data = 0, dim = c(nrow(geno), length(traits), length(models)),
                          dimnames = list(rownames(geno), traits, models))
  
  # User message
  cat(paste("\nMarker effects will be predicted using the following models:\n",paste(models, collapse="\n"),"\n",sep=""))
  
  
  # Create a temporary directory - delete this later
  dir.create(path = "tmp")
  
  ## Run predictions
  
  # Iterate over traits
  for (trt in traits) {
    
    # Print
    cat(sub("X", trt, "\nPredicting marker effects for trait X..."))
    
    # Format y, X
    mf <- model.frame(formula = reformulate("geno_id", trt), data = pheno)
    y <- model.response(mf)
    X <- model.matrix(~ 1, mf)
    M <- t(geno[, mf$geno_id, drop = FALSE])
    
    # First rrBLUP, if called for
    if ("rrBLUP" %in% models) {
      ans <- mixed.solve(y = y, Z = M, X = X, method = "REML")
      
      # Store
      marker_effects[names(ans$u), trt, "rrBLUP"] <- unname(ans$u)
      
    }
    
    # Now iterate over the bayesian models
    bayesian_models <- setdiff(models, "rrBLUP")
    
    if (length(bayesian_models) > 0) {
      
      for (mod in bayesian_models) {
        
        # Fit
        ETA <- list(
          list(X = X, model = "FIXED"),
          list(X = M, model = mod)
        )
        ans <- BGLR(y = y, response_type = "gaussian", ETA = ETA, nIter = params$nIter, burnIn = params$burnIn, 
                    saveAt = "tmp/temp", verbose = FALSE)
        
        marker_effects[names(ans$ETA[[2]]$b), trt, mod] <- unname(ans$ETA[[2]]$b)
        
      }
      
    }
    
    cat("done.")
    
  } # Stop trait loop
    
  # Remove the temporary file directory
  unlink(x = "tmp/", recursive = TRUE, force = TRUE)
  
  # Create the object
  output <- new(Class = "PopVar.me", ploidy = data@ploidy, map = data@map, geno.mat = data@geno.mat, haplo.mat = data@haplo.mat, 
                phased = data@phased, missing = data@missing, assume.inbred = data@assume.inbred, pheno = data@pheno, 
                fixed = data@fixed, traits = data@traits, marker.effects = marker_effects)
  
  return(output)

  
}
