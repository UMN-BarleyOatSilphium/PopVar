#' Predict genetic variance and genetic correlations in bi-parental populations deterministic equations
#' 
#' @param data A \code{PopVar.me} object from \code{\link{predict_marker_effects}}.
#' @param crossing.table A data frame of specific crosses to predict. Each row is a cross and the first two columns are the parents of the cross.
#' @param parents A character vector of potential parents. If passed, \code{crossing.table} is ignored and all pairwise crosses using \code{parents} will be predicted.
#' @param tail.p The percentile of the simulated progeny to be included into the calculation of \eqn{\mu}\emph{_sp} and correlated response. Default is \code{0.10}.
#' @param self.gen The number of selfing generations in the potential cross. Can be an integer or \code{Inf} for
#' recombinant inbreds. Note: \code{self.gen = 0} corresponds to an F1 population for outbreds and \code{self.gen = 1} corresponds to an F2 population.
#' @param DH Indicator if doubled-haploids are to be induced after the number of selfing generations indicated by
#' \code{self.gen}. For example, if \code{self.gen = 0} and \code{DH = TRUE}, then doubled-haploids are assumed
#' to be induced using gametes from F1 plants.
#' 
#' @details
#'
#' Predictions are based on the deterministic equations specified by Zhong and Jannink (2007), Allier et al. (2019),
#' and Neyhart et al. (2019).
#' 
#' 
#' @return
#' 
#' 
#' @examples 
#' 
#' \donttest{
#' 
#' # Read in phased cranberry marker data and phenotypes
#' geno <- cranberry_geno_phased
#' pheno <- cranberry_pheno
#' dat <- read_PopVar(geno = geno, pheno = pheno, n.traits = 6, phased = TRUE, 
#'                    assume.inbred = FALSE, dominance = TRUE, sep = "_")
#' 
#' # Compute marker effects using the rrBLUP model
#' data_me <- predict_marker_effects(data = dat, models = c("rrBLUP"))
#' 
#' # Simulate a crossing table
#' parents <- colnames(slot(data_me, "geno.mat"))
#' parent1 <- sample(parents, 10)
#' parent2 <- sample(setdiff(parents, parent1), 10)
#' cross_tab <- data.frame(parent1 = parent1, parent2 = parent2)
#' 
#' # Use example data to make predictions
#' out <- pop.predict2(G.in = G.in_ex_imputed, y.in = y.in_ex, map.in = map.in_ex, 
#'                     crossing.table = cross.tab_ex)
#'                     
#' # Provide a vector of parents to predict all possible crosses (some parents
#' # have missing phenotypic data)
#' out <- pop.predict2(G.in = G.in_ex_imputed, y.in = y.in_ex, map.in = map.in_ex, 
#'                     parents = y.in_ex$Entry[1:5])
#'                     
#' # Make predictions for 5 crosses with various levels of inbreeding
#' out_list <- lapply(X = 1:10, FUN = function(self.gen) {
#'   out <- pop.predict2(G.in = G.in_ex_imputed, y.in = y.in_ex, map.in = map.in_ex, 
#'                       crossing.table = cross.tab_ex[1:5,], self.gen = self.gen)
#'   out$self.gen <- self.gen
#'   out })
#'                
#' # Plot predictions of grain yield genetic variance over levels of inbreeding
#' dat <- do.call("rbind", lapply(out_list, subset, trait == "Yield"))
#' plot(pred_varG ~ self.gen, data = dat, type = "b", 
#'      subset = parent1 == parent1[1] & parent2 == parent2[1])
#' 
#' }
#'
#' @references
#'
#' Zhong, S., and J.-L. Jannink, 2007 Using quantitative trait loci results to discriminate among crosses on the basis of their
#' progeny mean and variance. Genetics 177: 567–576. https://doi.org/10.1534/ genetics.107.075358
#'
#' Allier, A., L. Moreau, A. Charcosset, S. Teyssèdre, and C. Lehermeier, 2019 Usefulness Criterion and Post-selection Parental
#' Contributions in Multi-parental Crosses: Application to Polygenic Trait Introgression. G3 9: 1469–1479.
#' doi: 10.1534/g3.119.400129
#'
#' Neyhart, J.L., A.J. Lorenz, and K.P. Smith, 2019 Multi-trait Improvement by Predicting Genetic Correlations in Breeding
#' Crosses. G3 9: 3153-3165. doi: 10.1534/g3.119.400406
#'
#' @importFrom qtl mf.h
#'
#' @export
#'
#'
predict_pop <- function(data, crossing.table, parents, tail.p = 0.1, self.gen = 0, DH = FALSE) {
  
  # Error handling
  stopifnot(inherits(data, "PopVar.me"))
  
  # Self gen cannot be negative
  stopifnot(self.gen >= 0)
  # DH must be logical
  stopifnot(is.logical(DH))

  # If the crossing table is not missing, check that the parents are in the G.in input
  if (!missing(crossing.table)) {
    parents <- unique(unlist(crossing.table))
    
    # Make sure the parent names are not factors
    crossing.table <- as.data.frame(sapply(X = crossing.table, as.character), stringsAsFactors = FALSE)
    
  } else {
    if (missing(parents)) stop("If no crossing.table is provided, a list of parents must be supplied.")
    
    parents <- sort(parents)
    # Create a crossing table with all possible parent combinations
    crossing.table <- as.data.frame(t(combn(x = parents, m = 2)), stringsAsFactors = FALSE)
    names(crossing.table) <- c("parent1", "parent2")
    
  }
  
  # Make sure the parents are genotyped
  haplo_mat <- data@haplo.mat
  genotyped_indiv <- colnames(data@geno.mat)
  if (any(!parents %in% genotyped_indiv)) stop("Individuals in 'parents' or 'crossing.table' must be genotyped.")
  
  
  ## Compute marker recombination frequencies
  # Split markers by chromosome
  gen_map <- data@map
  map_list <- split(gen_map, gen_map$chrom)
  markers_chr <- lapply(map_list, "[[", 1)
  # Calculate separate centimorgan distance matrices per chromosome
  chr_cM <- lapply(X = map_list, FUN = function(x) as.matrix(dist(x$cM)))
  # Convert to recombination distance
  chr_c <- lapply(X = chr_cM, FUN = mf.h)
  
  
  
  # Calculate the LD covariance coefficient
  # If self.gen == 0, DH does not matter; it's the same coefficient
  if (self.gen == 0) {
    chr_covar <- lapply(X = chr_c, FUN = function(c) 1 - (2 * c))
  
  } else if (DH & self.gen > 0) {
    covar_selfing <- lapply(X = chr_c, FUN = function(c) {
      base <- 0.5 * (1 - (2 * c))
      Reduce(f = `+`, x = lapply(X = seq(self.gen), FUN = function(k) base ^ k)) })
    
    covar_final <- lapply(X = chr_c, FUN = function(c) (0.5 * (1 - (2 * c))) ^ self.gen)
    
    # Sum
    chr_covar <- mapply(FUN = `+`, covar_selfing, covar_final)
    
  } else if (!DH & is.finite(self.gen)) {
    chr_covar <- lapply(X = chr_c, FUN = function(c) {
      base <- 0.5 * (1 - (2 * c))
      Reduce(f = `+`, x = lapply(X = seq(self.gen), FUN = function(k) base ^ k)) })
    
  } else if (!DH & is.infinite(self.gen)) {
    chr_covar <- lapply(X = chr_c, FUN = function(c) (1 - (2 * c)) / (1 + (2 * c)))
    
  }
  
  
  # Predict breeding values and total genotypic values depending on whether dominance was specified
  if (data@dominance) {
    
  }
    
    
    
  
  
  
  
  ## Predicted genotypic value of all genotypes + grand mean
  pgvs <- (M %*% mar_eff_mat) + matrix(mar_beta_mat, ncol = ncol(mar_eff_mat), nrow = nrow(M), byrow = TRUE)
  
  
  ## Calculate the pairwise product of marker effects for each trait, separated by chromosome
  ## Then multiply by the LD covariance
  intra_trait_covar <- list()
  
  # Iterate over traits
  for (i in seq(ncol(mar_eff_mat))) {
    
    # Split the marker effect matrix by chromosome and calculate pairwise product
    mar_eff_prod_chr <- lapply(X = markers_chr, function(marks) tcrossprod(mar_eff_mat[marks, i, drop = F]))
    
    # The covariance is the QTL effect product multiplied by the expected D
    intra_trait_covar[[trait_names[i]]] <- mapply(mar_eff_prod_chr, chr_covar, FUN = `*`, SIMPLIFY = FALSE)
    
  }
  
  
  if (n_traits > 1) {
    
    ## Calculate the pairwise product of marker effects for each pair of traits
    # First create a vector of trait index combinations
    trait_ind_combn <- t(combn(x = seq(n_traits), m = 2, simplify = TRUE))
    trait_combn_name <- combn(x = trait_names, m = 2, paste0, collapse = ":")
    
    # Empty matrix for correlations
    trait_corG_mat <- matrix(data = NA, nrow = n_traits, ncol = n_traits, dimnames = list(trait_names, paste0("cor_w_", trait_names)))
    # Convert to a distance version for easy addition
    trait_corG_dist <- as.dist(trait_corG_mat)
    
    # List of inter-trait covariances
    inter_trait_covar <- list()
    
    # Iterate over trait indices combinations
    for (i in  seq(nrow(trait_ind_combn))) {
      
      # Get the marker effects for that combination
      mar_eff_combn <- mar_eff_mat[,trait_ind_combn[i,], drop = FALSE]
      
      # Calculate the pairwise product by chromosome
      mar_eff_prod_chr <- lapply(X = markers_chr, function(marks) tcrossprod(mar_eff_combn[marks,1,drop = FALSE], mar_eff_combn[marks,2,drop = FALSE]))
      
      # The covariance is the QTL effect product multiplied by the expected D
      inter_trait_covar[[trait_combn_name[i]]] <- mapply(mar_eff_prod_chr, chr_covar, FUN = `*`, SIMPLIFY = FALSE)
      
    }
    
  } else {
    inter_trait_covar <- NULL
    
  }
  
  
  # Determine the k_sp from the the tail.p
  k_sp <- dnorm(x = qnorm(p = 1 - tail.p)) / tail.p
  
  
  ## Verify that all parents are homozygous for all markers
  parents_M <- M1[parents,,drop = FALSE]
  
  # Make sure the markers are coded correctly
  if (!all(parents_M %in% c(-1, 1)))
    stop("This method assumes fully inbred parents. Therefore, marker genotypes other than -1 and 1 are not supported. Please
    edit the marker data and try again.")
  
  ## Create a list to store dfs
  cross_predictions <- vector("list", length = nrow(crossing.table))
  
  ## Iterate over the pairs of parents
  for (j in seq(nrow(crossing.table))) {
    
    # Character vector of the two parents
    pars <- as.character(crossing.table[j,1:2])
    # Cross mean prediction
    pred_mu_j <- colMeans(pgvs[pars,,drop = FALSE])
    # Genetic variance - use the pred_mu_j vector as a template
    pred_varG_j <- pred_mu_j
    
    ## Subset the genotype matrix using the parents
    par_geno <- M1[pars,,drop = FALSE]
    
    # Which markers are segregating?
    mar_seg <- names(which(colMeans(par_geno) == 0))
    # If no markers are segregating, the variance is 0
    if (length(mar_seg) == 0) {
      pred_varG_j[] <- 0
      pred_corG_mat <- if (n_traits > 1) trait_corG_mat else NULL
      pred_cor_musp_low <- pred_cor_musp_high <- if (n_traits > 1) pred_mu_j else NULL
      
    } else {
      
      # Split by chromsome
      mar_seg_chr <- lapply(markers_chr, intersect, mar_seg)
      
      # Find the parent 1 genotype of those markers and take the crossproduct for multiplication
      par1_mar_seg <- crossprod(par_geno[1,mar_seg, drop = FALSE])
      # Split by chromosome
      par1_mar_seg_chr <- lapply(mar_seg_chr, function(marks) par1_mar_seg[marks, marks, drop = FALSE])
      
      
      ## Predictions
      
      # Iterate over traits
      # Note that the QTL covariance matrix includes the variance of each QTL on the diagonal, so the sum of the matrix
      # is the variance + 2 * covariance
      for (i in seq(n_traits)) {
        pred_varG_j[i] <- sum(mapply(par1_mar_seg_chr, intra_trait_covar[[i]], FUN = function(x, y)
          sum(x * y[colnames(x), colnames(x)])))
      }
      
      
      # Genetic correlations between traits, if more than one trait
      if (!is.null(inter_trait_covar)) {
        
        ## Iterate over trait combinations
        for (i in seq_along(inter_trait_covar)) {
          
          ## Calculate the covariance between the pair of traits
          trait_pair_cov <- sum(mapply(par1_mar_seg_chr, inter_trait_covar[[i]], FUN = function(x, y) sum(x * y[colnames(x), colnames(x)])))
          
          # Subset predicted variance
          trait_pair_varG <- pred_varG_j[trait_ind_combn[i,]]
          
          # Calculate correlation and save
          trait_corG_dist[i] <- trait_pair_cov / prod(sqrt(trait_pair_varG))
          
        }
        
        ## Convert distance object to matrix
        pred_corG_mat <- as.matrix(trait_corG_dist)
        dimnames(pred_corG_mat) <- dimnames(trait_corG_mat)
        diag(pred_corG_mat) <- NA
        
        ## Calculate correlated progeny mean
        response_trait_varG <- matrix(pred_varG_j, nrow = length(pred_varG_j), ncol = length(pred_varG_j), byrow = TRUE)
        correlated_response <- k_sp * pred_corG_mat * sqrt(response_trait_varG)
        pred_mu_j_mat <- matrix(pred_mu_j, nrow = length(pred_mu_j), ncol = length(pred_mu_j), byrow = TRUE)
        pred_cor_musp_low <- pred_mu_j_mat - correlated_response
        pred_cor_musp_high <- pred_mu_j_mat + correlated_response
        
        # Change names
        colnames(pred_cor_musp_low) <- paste0("pred_cor_musp_low_", trait_names)
        colnames(pred_cor_musp_high) <- paste0("pred_cor_musp_high_", trait_names)
        
      } else {
        pred_corG_mat <- pred_cor_musp_low <- pred_cor_musp_high <- NULL
        
      }
      
    }
    
    
    ## Save the results as a data.frame
    cross_predictions[[j]] <- data.frame(parent1 = pars[1], parent2 = pars[2], trait = trait_names,
                                         cbind(pred_mu = pred_mu_j, pred_varG = pred_varG_j, pred_corG_mat, pred_cor_musp_low, 
                                               pred_cor_musp_high),
                                         stringsAsFactors = FALSE, row.names = NULL)
    
  } # End loop
  
  ## Bind rows
  cross_predictions1 <- do.call("rbind", cross_predictions)
  
  ## Calculate response predictions (superior progeny, correlated response, etc.)
  pred_response <- (k_sp * sqrt(cross_predictions1$pred_varG))
  
  # Superior progeny mean
  cross_predictions1[["pred_musp_low"]] <- cross_predictions1$pred_mu - pred_response
  cross_predictions1[["pred_musp_high"]] <- cross_predictions1$pred_mu + pred_response
  
  ## Re-order columns
  cor_W_cols <- grep(pattern = paste0(trait_names, collapse = "|"), x = names(cross_predictions1))
  cross_predictions2 <- cross_predictions1[, c(setdiff(seq(ncol(cross_predictions1)), cor_W_cols), cor_W_cols)]
  
  ## Return the predictions
  return(cross_predictions2)
  
} # Close function


  
  
  
  
  
  
  