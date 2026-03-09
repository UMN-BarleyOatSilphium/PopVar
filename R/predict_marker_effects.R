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
#' geno <- cranberry_geno_phased
#' pheno <- cranberry_pheno
#' dat <- read_PopVar(geno = geno, pheno = pheno, n.traits = 6, phased = TRUE, 
#'                    assume.inbred = FALSE, dominance = TRUE, sep = "_")
#' 
#' # Compute marker effects using the rrBLUP and RKHS models
#' data_me <- predict_marker_effects(data = dat, models = c("rrBLUP", "BRR"))
#' 
#' @importFrom rrBLUP mixed.solve
#' @import sommer
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
  add_marker_effects <- array(data = 0, dim = c(nrow(geno), length(traits), length(models)),
                              dimnames = list(rownames(geno), traits, models))
  allele_sub_effects <- dom_deviation_effects <- dom_marker_effects <- add_marker_effects
  
  # User message
  cat(paste("\nMarker effects will be predicted using the following models:\n",paste(models, collapse="\n"),"\n",sep=""))
  
  
  # Create a temporary directory - delete this later
  dir.create(path = "tmp", showWarnings = FALSE)
  
  ## Run predictions
  
  # Iterate over traits
  for (trt in traits) {
    
    # Print
    cat(sub("X", trt, "\nPredicting marker effects for trait X..."))
    
    # Format y, X
    mf <- model.frame(formula = reformulate("geno_id", trt), data = pheno)
    y <- model.response(mf)
    X <- model.matrix(~ 1, mf)
    # M for mmes and BGLR must be -1, 0, 1
    # M <- t(geno[, mf$geno_id, drop = FALSE] - 1)
    # Or should it be 0, 1, 2?
    M <- t(geno[, mf$geno_id, drop = FALSE])
    
    # Center the matrices
    p_j <- data@allele.freq
    Z <- M - matrix(2 * p_j, nrow = nrow(M), ncol = ncol(M), byrow = TRUE)
    M1 <- (M - 1) - matrix(2 * p_j, nrow = nrow(M), ncol = ncol(M), byrow = TRUE)
    
    if (data@dominance) {
      # Create a dominance matrix
      R <- geno[, mf$geno_id, drop = FALSE]
      R[R != 1] <- 0
      R <- t(R)
      
      q_j <- 1 - p_j
      R1 <- R - matrix(2 * p_j * q_j, nrow = nrow(R), ncol = ncol(R), byrow = TRUE)
      
    }
      
    
    # First rrBLUP, if called for
    if ("rrBLUP" %in% models) {
      # Single variance component if no dominance
      if (!data@dominance) {
        ans <- mixed.solve(y = y, Z = Z, X = X, method = "REML")
        # ans1 <- mixed.solve(y = y, Z = M, X = X, method = "REML")
        # Store
        add_marker_effects[names(ans$u), trt, "rrBLUP"] <- unname(ans$u)
        allele_sub_effects[names(ans$u), trt, "rrBLUP"] <- unname(ans$u)
        
      } else {
        mf1 <- mf
        # Add directional dominance effect
        mf1$f <- rowMeans(M != 0)
        X1 <- model.matrix(~ 1 + f, mf1)
        
        
        # #### TEST ###
        # Ztest <- Z[,data@map$marker[data@map$chrom == "chr01"]]
        # R1test <- R1[,data@map$marker[data@map$chrom == "chr01"]]
        # 
        # # SOMMER for direct marker effects
        # ans_a <- mmes(fixed = reformulate(termlabels = c("1", "f"), response = trt),
        #               random = ~ vsm(ism(Ztest)) + vsm(ism(R1test)), data = mf1, verbose = FALSE, dateWarning = FALSE)
        #             
        # 
        # # EMMREML
        # ans_b <- emmremlMultiKernel(y = y, X = X1, Zlist = list(a = Ztest, d = R1test), Klist = list(a = diag(ncol(Ztest)), d = diag(ncol(Ztest))))
        # 
        # # Compare
        # plot(ans_a$uList[[1]], ans_b$uhat[seq(1, ncol(Ztest))]); abline(a = 0, b = 1)
        # plot(ans_a$uList[[2]], ans_b$uhat[seq(ncol(Ztest) + 1, (ncol(Ztest) * 2))]); abline(a = 0, b = 1)
        # 
        # 
        # # Alternative parameterization
        # mf1$geno_id <- as.character(mf1$geno_id)
        # mf1$geno_id <- factor(mf1$geno_id, levels = mf1$geno_id)
        # mf1$geno_id1 <- mf1$geno_id
        # 
        # 
        # MMT <- tcrossprod(Ztest)
        # MMTinv <- solve(MMT)
        # MMTinv <- t(Ztest) %*% MMTinv
        # 
        # RRT <- tcrossprod(R1test)
        # RRTinv <- solve(RRT)
        # RRTinv <- t(R1test) %*% RRTinv
        # 
        # ans1 <- mmes(fixed = reformulate(termlabels = c("1", "f"), response = trt),
        #              random = ~ vsm(ism(geno_id), Gu = MMT) + vsm(ism(geno_id1), Gu = RRT),
        #              data = mf1, verbose = FALSE, dateWarning = FALSE)
        # 
        # # Compare
        # add.me.part <- MMTinv %*% as.matrix(ans1$uList$`vsm(ism(geno_id), Gu = MMT`)
        # dom.me.part <- RRTinv %*% as.matrix(ans1$uList$`vsm(ism(geno_id1), Gu = RRT`)
        # 
        # plot(ans_a$uList[[1]], add.me.part); abline(a = 0, b = 1)
        # plot(ans_a$uList[[2]], dom.me.part); abline(a = 0, b = 1)
        
        ##### End test ### 
        
        # Alternative parameterization
        mf1$geno_id <- as.character(mf1$geno_id)
        mf1$geno_id <- factor(mf1$geno_id, levels = mf1$geno_id)
        mf1$geno_id1 <- mf1$geno_id
        
        MMT <- tcrossprod(Z)
        MMTinv <- solve(MMT)
        MMTinv <- t(Z) %*% MMTinv
        
        RRT <- tcrossprod(R1)
        RRTinv <- solve(RRT)
        RRTinv <- t(R1) %*% RRTinv

        try_ans <- try({
          ans1 <- mmes(fixed = reformulate(termlabels = c("1", "f"), response = trt),
                       random = ~ vsm(ism(geno_id), Gu = MMT) + vsm(ism(geno_id1), Gu = RRT),
                       data = mf1, verbose = FALSE, dateWarning = FALSE)
        }, silent = TRUE)
        
        # If error - fit only additive
        if (inherits(try_ans, "try-error")) {
          ans1 <- mmes(fixed = reformulate(termlabels = c("1", "f"), response = trt),
                       random = ~ vsm(ism(geno_id), Gu = MMT),
                       data = mf1, verbose = FALSE, dateWarning = FALSE)
          
          add.me.part <- MMTinv %*% as.matrix(ans1$uList$`vsm(ism(geno_id), Gu = MMT`)
          dom.me.part <- add.me.part
          dom.me.part[,] <- 0
          
        } else {
          add.me.part <- MMTinv %*% as.matrix(ans1$uList$`vsm(ism(geno_id), Gu = MMT`)
          dom.me.part <- RRTinv %*% as.matrix(ans1$uList$`vsm(ism(geno_id1), Gu = RRT`)
          
        }

        # Store
        add_marker_effects[row.names(add.me.part), trt, "rrBLUP"] <- as.vector(add.me.part)
        dom_marker_effects[row.names(dom.me.part), trt, "rrBLUP"] <- as.vector(dom.me.part)
        
        ## Convert marker effects a and d to substitution effects and dominance deviations
        # Get the regression coefficient on percent homozygosity
        f_coef <- ans1$b[2,1]
        # Convert dom effects to dom deviations
        dom_dev <- as.vector(dom.me.part) - (f_coef / length(dom.me.part))
        # Compute substitution effects
        allele_sub <- as.vector(add.me.part) + (dom_dev * (q_j - p_j))
        
        # Add these to the matrices
        allele_sub_effects[row.names(add.me.part), trt, "rrBLUP"] <- unname(allele_sub)
        dom_deviation_effects[row.names(add.me.part), trt, "rrBLUP"] <- dom_dev
        
        # # Verify
        # R2 <- apply(X = M, MARGIN = 1, FUN = function(ind) {
        #   ind[ind == 1] <- 2 * (p_j[ind == 1] * q_j[ind == 1])
        #   ind[ind == 2] <- -2 * (q_j[ind == 2]^2)
        #   ind[ind == 0] <- -2 * (p_j[ind == 0]^2)
        #   ind
        # })
        # R2 <- t(R2)
        # RRT2 <- tcrossprod(R2)
        # RRT2inv <- solve(RRT2)
        # RRT2inv <- t(R2) %*% RRT2inv
        # ans2 <- mmes(fixed = reformulate(termlabels = c("1"), response = trt),
        #              random = ~ vsm(ism(geno_id), Gu = MMT) + vsm(ism(geno_id1), Gu = RRT2),
        #              data = mf1, verbose = FALSE, dateWarning = FALSE)
        # a_est <- MMTinv %*% as.matrix(ans2$uList$`vsm(ism(geno_id), Gu = MMT`)
        # # Compare allele effects from model 1 and model 2
        # plot(a_est[,1], add.me.part[,1]); abline(a = 0, b = 1)
        # # Compare allele effects from model 1 with recomputed allele substitution effects
        # plot(a_est[,1], allele_sub); abline(a = 0, b = 1)
        # 
        # plot(as.vector(add.me.part), unname(ans$u))
        
      }
      
 
      
    }
    
    # Now iterate over the bayesian models
    bayesian_models <- setdiff(models, "rrBLUP")
    
    if (length(bayesian_models) > 0) {
      
      for (mod in bayesian_models) {
        
        # Additive or dominance
        if (!data@dominance) {
        
          # Fit
          ETA <- list(
            list(X = X, model = "FIXED"),
            list(X = Z, model = mod)
          )
          ans <- BGLR(y = y, response_type = "gaussian", ETA = ETA, nIter = params$nIter, burnIn = params$burnIn, 
                      saveAt = "tmp/temp", verbose = FALSE)
          
          add_marker_effects[names(ans$ETA[[2]]$b), trt, mod] <- unname(ans$ETA[[2]]$b)
          allele_sub_effects[names(ans$ETA[[2]]$b), trt, mod] <- unname(ans$ETA[[2]]$b)
          
          
        } else {
          
          X1 <- model.matrix(~ 1 + f, mf1)
          
          # Fit
          ETA <- list(
            fixed = list(X = X1, model = "FIXED"),
            add = list(X = Z, model = mod),
            dom = list(X = R1, model = mod)
          )
          ans <- BGLR(y = y, response_type = "gaussian", ETA = ETA, nIter = params$nIter, burnIn = params$burnIn, 
                      saveAt = "tmp/temp", verbose = FALSE)
          
          add_marker_effects[names(ans$ETA[[2]]$b), trt, mod] <- unname(ans$ETA[[2]]$b)
          dom_marker_effects[names(ans$ETA[[3]]$b), trt, mod] <- unname(ans$ETA[[3]]$b)
          
          
          ## Convert marker effects a and d to substitution effects and dominance deviations
          # Get the regression coefficient on percent homozygosity
          f_coef <- ans$ETA$fixed$b[2]
          p <- length(unname(ans$ETA[[3]]$b))
          # Convert dom effects to dom deviations
          dom_dev <- unname(ans$ETA[[3]]$b) - (f_coef / p)
          # Compute substitution effects
          allele_sub <- unname(ans$ETA[[2]]$b) + (dom_dev * (q_j - p_j))
          
          # Add these to the matrices
          allele_sub_effects[names(ans$ETA[[2]]$b), trt, mod]  <- unname(allele_sub)
          dom_deviation_effects[names(ans$ETA[[3]]$b), trt, mod]  <- dom_dev
          
        }
        
      }
      
    }
    
    cat("done.")
    
  } # Stop trait loop
    
  # Remove the temporary file directory
  unlink(x = "tmp/", recursive = TRUE, force = TRUE)

  marker_effects <- list(add.marker.effects = add_marker_effects, dom.marker.effects = dom_marker_effects, 
                         allele.sub.effects = allele_sub_effects, dom.dev.effects = dom_deviation_effects)
  
  output <- new(Class = "PopVar.me", ploidy = data@ploidy, map = data@map, geno.mat = data@geno.mat,  
                haplo.mat = data@haplo.mat, allele.freq = data@allele.freq, phased = data@phased, missing = data@missing, 
                assume.inbred = data@assume.inbred, dominance = data@dominance, pheno = data@pheno, 
                fixed = data@fixed, traits = data@traits, marker.effects = marker_effects)
  
  return(output)

  
}
