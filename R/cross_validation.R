#' Estimate genomewide prediction accuracy using cross-validation
#' 
#' 
#' @param data A \code{PopVar.data} object from \code{read_PopVar}.
#' @param group The name of a column in the phenotype data for grouped cross-validation (e.g. environment, location, etc.). For example to run leave-one-environment-out cross-validation, set \code{method = "loo"} and \code{group = "env"}.
#' @param models The list of models to run. See Details for supported models.
#' @param method The cross-validation method. Either fractional (\code{"frac"}), k-fold (\code{"kfold"}), or leave-one-out \code{"loo"}.
#' @param p.test For fractional cross-validation, the proportion of individuals to designate as the test set per replicate. \code{1 - p.test} is the proportion of individuals used for training.
#' @param k For k-fold cross-validation, the number of folds.
#' @param n.reps The number of complete replicates of fractional or k-fold cross-validation.
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
#' @return A list with a summary table and plot of cross-validation results.
#' 
#' @examples
#' # Read in unphased cranberry marker data and phenotypes
#' geno <- cranberry_geno_phased
#' pheno <- cranberry_pheno
#' dat <- read_PopVar(geno = geno, pheno = pheno, n.traits = 6, phased = TRUE, 
#'                    assume.inbred = FALSE, dominance = TRUE, sep = "_")
#' 
#' # Run k-fold cross-validation using the rrBLUP and RKHS models
#' cv_res <- cross_validate(data = dat, models = c("rrBLUP", "BRR"), method = "kfold", n.reps = 10)
#' 
#' @importFrom rrBLUP mixed.solve
#' @import sommer
#' @importFrom BGLR BGLR
#' @import ggplot2
#' 
#' 
#' @export
#' 
cross_validate <- function(data, group = NULL, models = c("rrBLUP", "BayesA", "BayesB", "BayesC", "BL", "BRR"),
                           method = c("frac", "kfold", "loo"), params = list(nIter = 500, burnIn = 50), 
                           p.test = 0.3, k = 5, n.reps = 25) {
  
  # Error handling
  stopifnot(inherits(data, "PopVar.data"))
  models <- match.arg(models, several.ok = TRUE)
  stopifnot(is.list(params))
  required_param_names <- c("nIter", "burnIn")
  if (!all(required_param_names %in% names(params))) stop("The following must be in params: ", paste0(required_param_names, collapse = ", "))
  
  stopifnot(is.numeric(p.test))
  stopifnot(is.numeric(k))
  stopifnot(is.numeric(n.reps))
  
  stopifnot(p.test > 0 & p.test < 1)
  stopifnot(k > 1)
  stopifnot(n.reps > 0)
  
  # Integerize
  k <- as.integer(k)
  n.reps <- as.integer(n.reps)
  
  # Get the genotype data
  geno <- data@geno.mat
  # Get phenotypes
  pheno <- data@pheno
  pheno$geno_id <- factor(pheno$geno_id, levels = colnames(geno))
  # Get traits
  traits <- data@traits
  
  # Create additive and (if flagged) dominance relationship matrices
  M <- t(geno - 1)
  Amat <- A.mat(X = M, min.MAF = 0, return.imputed = FALSE)
  
  # Center the matrices
  p_j <- data@allele.freq
  M1 <- M - matrix(2 * p_j, nrow = nrow(M), ncol = ncol(M), byrow = TRUE)
  
  if (data@dominance) {
    Dmat <- D.mat(X = M, nishio = FALSE, min.MAF = 0, return.imputed = FALSE)
    R <- M + 1
    R[R != 1] <- 0
    est.f <- rowMeans(M == 1)
    pheno$f <- est.f[match(pheno$geno_id, names(est.f))]

    q_j <- 1 - p_j
    R1 <- R - matrix(2 * p_j * q_j, nrow = nrow(R), ncol = ncol(R), byrow = TRUE)
    
  }
  
  # User message
  method_print <- if (method == "kfold") "k-fold" else if (method == "frac") "fractional" else "leave-one-out"
  message <- sub("X", method_print, "\nRunning Y reps of X cross-validation using the following models:\n")
  message <- sub("Y", n.reps, message)
  cat(paste(message, paste(models, collapse = "\n"), "\n", sep = ""))
  
  # Create a temporary directory - delete this later
  dir.create(path = "tmp", showWarnings = FALSE)
  
  
  # Create randomizations
  cv_rand <- list()
  
  # Iterate over traits
  for (trt in traits) {
    # Which are na
    idx_na <- which(!is.na(pheno[[trt]]))
    pheno_idx <- pheno[idx_na, c(colnames(pheno)[1], trt)]
    
    # Randomize
    if (method == "frac") {
      train_test_idx_list <- replicate(n = n.reps, {
        test <- sort(sample(x = idx_na, size = ceiling(p.test * length(idx_na))))
        train <- setdiff(idx_na, test)
        list(train = train, test = test)
      }, simplify = FALSE)
      
    } else if (method == "kfold") {
      train_test_idx_list <- replicate(n = n.reps, {
        # Break up idx_na into folds
        idx_na_folds <- split(x = idx_na, f = sample(cut(x = idx_na, breaks = k)))
        # Iterate and create training/testing sets
        fold_train_test <- lapply(X = idx_na_folds, FUN = function(idxx) {
          train <- setdiff(idx_na, idxx)
          test <- idxx
          list(train = train, test = test)
        })
        unname(fold_train_test)
      }, simplify = FALSE)
      
    } else {
      train_test_idx_list <- list(lapply(X = idx_na, FUN = function(idxx) {
        train <- setdiff(idx_na, idxx)
        test <- idxx
        list(train = train, test = test)
      }))
      
    }
    
    # Add to the list
    cv_rand[[trt]] <- train_test_idx_list
    
  }
  
  
  ## Run cross-validation
  
  # A list to store cv results
  cv_out <- list()
  
  # Iterate over the randomization
  for (i in seq_along(cv_rand)) {
    
    trt <- names(cv_rand)[i]
    rand_trt <- cv_rand[[i]]

    # Print
    cat(sub("X", trt, "\nRunning cross-validation for trait X..."))
  
    # List to store accuracy results
    acc_list <- list()
      
    # Reformat the rand_trt if the method is frac
    if (method == "frac") {
      rand_trt <- lapply(rand_trt, list)
    }
    
    # Iterate over reps
    for (ii in seq_along(rand_trt)) {
      rand_trt_ii <- rand_trt[[ii]]
      
      # List to store fold prediction results
      fold_pred_list <- list()
      
      # Iterate over folds
      for (f in seq_along(rand_trt_ii)) {
        rand_f <- rand_trt_ii[[f]]
        
        # A matrix to store predictions
        preds_f <- matrix(data = 0, nrow = length(rand_f$test), ncol = length(models) + 1, dimnames = list(pheno$geno_id[rand_f$test], c("values", models)))
        
        # Format y, X
        test <- pheno[rand_f$test, trt, drop = FALSE]
        row.names(test) <- pheno$geno_id[rand_f$test]
        preds_f[,"values"] <- test[[trt]]
        
        form <- if (data@dominance) reformulate(c("geno_id", "f"), trt) else reformulate("geno_id", trt)
        mf <- model.frame(formula = form, data = pheno[rand_f$train, ])
        y <- model.response(mf)
        X <- model.matrix(~ 1, mf)
        Z <- model.matrix(~ -1 + geno_id, mf)
        
        # First rrBLUP, if called for
        if ("rrBLUP" %in% models) {
          
          # Single variance component if no dominance
          if (!data@dominance) {
            ans <- mixed.solve(y = y, Z = Z, X = X, K = Amat, method = "REML")
            blups <- ans$u
            
          } else {
            mf1 <- mf
            mf1$geno_id1 <- mf1$geno_id

            try_ans <- try({
              ans1 <- mmes(fixed = reformulate(termlabels = c("1"), response = trt),
                           random = ~ vsm(ism(geno_id), Gu = Amat) + vsm(ism(geno_id1), Gu = Dmat),
                           data = mf1, verbose = FALSE, dateWarning = FALSE)
            }, silent = TRUE)
            
            blups <- ans1$uList[[1]] + ans1$uList[[2]]
            blups <- blups[,1]

          }
          
          # Add predictions
          preds_f[,"rrBLUP"] <- blups[row.names(preds_f)]
          
        }
        
        # Now iterate over the bayesian models
        bayesian_models <- setdiff(models, "rrBLUP")
        
        if (length(bayesian_models) > 0) {
          
          for (mod in bayesian_models) {
            
            # Additive or dominance
            if (!data@dominance) {
              
              # Fit
              ETA <- list(
                fixed = list(X = X, model = "FIXED"),
                additive = list(X = M1[mf$geno_id, , drop = FALSE], model = mod)
              )
              ans <- BGLR(y = y, response_type = "gaussian", ETA = ETA, nIter = params$nIter, burnIn = params$burnIn, 
                          saveAt = "tmp/temp", verbose = FALSE)
              
              # Get marker effects
              mar_eff1 <- ans$ETA$additive$b
              blups <- M1 %*% mar_eff1
              blups <- blups[,1]
              

            } else {
              
              # X1 <- model.matrix(~ 1 + f, mf1)
              X1 <- X
              
              # Fit
              ETA <- list(
                fixed = list(X = X1, model = "FIXED"),
                additive = list(X = M1[mf$geno_id, , drop = FALSE], model = mod),
                dominance = list(X = R1[mf$geno_id, , drop = FALSE], model = mod)
              )
              ans <- BGLR(y = y, response_type = "gaussian", ETA = ETA, nIter = params$nIter, burnIn = params$burnIn, 
                          saveAt = "tmp/temp", verbose = FALSE)
              
              # Get marker effects
              mar_eff1 <- ans$ETA$additive$b
              mar_eff2 <- ans$ETA$dominance$b
              blups <- (M1 %*% mar_eff1) + (R1 %*% mar_eff2)
              blups <- blups[,1]
              
            }
            
            # Add predictions
            preds_f[,mod] <- blups[row.names(preds_f)]
            
          }
          
        }
        
        # Add predictions to the fold list
        fold_pred_list[[f]] <- preds_f
        
      } # End fold loop
      
      # Bind rows
      fold_pred_mat <- do.call(rbind, fold_pred_list)
      acc <- cor(fold_pred_mat, use = "pairwise.complete.obs")
      
      
      # Compute accuracy; add to the list
      acc_list[[ii]] <- data.frame(rep = ii, stat = c("accuracy"), acc["values",models, drop = FALSE], row.names = NULL)
      
    } # End of rep iterator
      
    # Bind the accuracy results
    acc_df <- do.call(rbind, acc_list)
    acc_df <- cbind(trait = trt, acc_df)
    
    # Store results
    cv_out[[trt]] <- acc_df
    
    cat("done.")
    
  } # Stop trait loop
  
  # Combine the output list
  cv_out_df <- do.call(rbind, cv_out)
  row.names(cv_out_df) <- NULL
  
  # Create a plot
  cv_out_df1 <- reshape(data = cv_out_df, varying = models, v.names = "stat.value", timevar = "model", times = models, direction = "long")
  
  g <- ggplot(data = cv_out_df1, aes(x = trait, y = stat.value, fill = model)) +
    geom_boxplot(alpha = 0.75) +
    facet_grid(stat ~ ., scales = "free_y", switch = "y") +
    scale_y_continuous(name = NULL, breaks = pretty, guide = guide_axis(cap = TRUE)) +
    theme_classic() +
    theme(strip.placement = "outside", strip.background.y = element_blank())

  # Compute a summary
  cv_summ <- aggregate(stat.value ~ model + trait + stat, data = cv_out_df1, FUN = function(x) c(mean = mean(x), sd = sd(x)))
  cv_summ <- cbind(cv_summ[1:3], cv_summ$stat.value)
  
  # Remove the temporary file directory
  unlink(x = "tmp/", recursive = TRUE, force = TRUE)
  
  # Return the summaries
  output <- list(summary = cv_summ, plot = g)
  return(output)
  
}