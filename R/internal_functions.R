#' Internal functions
#' 
#' @name internal
#' 
#' @description 
#' Internal functions generally not to be called by the user.
#' 
#' @param crossing.table The crossing table.
#' @param par.entries The parent entries.
#' @param crossing.mat The crossing matrix.
#' @param GEBVs The genomic estimated breeding values.
#' @param tail.p The proportion from the population to select.
#' @param G The marker genotypes
#' @param y.CV The phenotypes for cross-validation.
#' @param G.CV The marker genotypes for cross-validation.
#' @param models.CV The models for cross-validation.
#' @param frac.train.CV The fraction of data to use as training data in cross-validation.
#' @param nCV.iter.CV The number of iterations of cross-validation.
#' @param burnIn.CV The burn-in number for cross-validation.
#' @param nIter.CV The number of iterations for Bayesian models in cross-validation.
#' @param nFold.CV The number of folds in k-fold cross-validation.
#' @param nFold.CV.reps The number of replications of k-fold cross-validation.
#' @param M The marker matrix.
#' @param y.df The phenotype data.
#' @param models The models to use.
#' @param nIter The number of iterations.
#' @param burnIn The burn-in rate.
#' 
#' @importFrom stats cor sd
#' @importFrom BGLR BGLR
#' 

#' 
#' @rdname internal
#' 
par_position <- function(crossing.table, par.entries) { # Used when a crossing table is defined
  
  par.pos <- matrix(nrow = nrow(crossing.table), ncol = 2)
  crosses.possible <- matrix(nrow = nrow(crossing.table), ncol = 2)
  for(i in 1:nrow(crossing.table)){
    par1 <- as.character(crossing.table[i,1])
    par2 <- as.character(crossing.table[i,2])
    
    if(par1 %in% par.entries & par2 %in% par.entries){
      par.pos[i,1] <- which(par.entries == par1)
      par.pos[i,2] <- which(par.entries == par2)
      crosses.possible[i,1] <- par1
      crosses.possible[i,2] <- par2  
    }
  }
  
  par.pos <- par.pos[which(!is.na(par.pos[,1])), ]
  crosses.possible <- crosses.possible[which(!is.na(crosses.possible[,1])), ]
  
  for.dup <- paste(par.pos[,1], ".", par.pos[,2], sep=""); duplicated <- which(duplicated(for.dup))
  if(length(duplicated) > 0){
    par.pos <- par.pos[-duplicated, ]
    crosses.possible <- crosses.possible[-duplicated, ]
  }
  
  return(list(parent.positions=par.pos, crosses.possible=crosses.possible))
}

#'
#' @rdname internal
#' 
par_name <- function(crossing.mat, par.entries){ ## Used when all combinations of parents are crossed
  crosses.possible <- matrix(nrow = nrow(crossing.mat), ncol = 2)
  for(i in 1:nrow(crossing.mat)){
    crosses.possible[i,1] <- par.entries[crossing.mat[i,1]]
    crosses.possible[i,2] <- par.entries[crossing.mat[i,2]]
  }
  return(crosses.possible)
}

#'
#' @rdname internal
#' 
tails <- function(GEBVs, tail.p){ #Calculates means of tails; set tail.p to the proportion of the populaiton you want to take the mean of, default is 10%
  u.top <- mean(GEBVs[which(GEBVs >= quantile(GEBVs, 1-tail.p))], na.rm=T)
  u.bot <- mean(GEBVs[which(GEBVs <= quantile(GEBVs, tail.p))], na.rm=T)
  
  return(rbind(u.top, u.bot))
}

#'
#' @rdname internal
#' 
maf_filt <- function(G){
  G_noNA <- sum(!is.na(G), na.rm = TRUE)
  freq1 <- sum(G == 1, na.rm = TRUE) / G_noNA #+ .5*(sum(G == 0, na.rm = TRUE) / G_noNA)
  min(freq1, 1 - freq1)
}


#test.map <- qtl::sim.map(len = c(50,50), n.mar = c(100,100), anchor.tel = FALSE, include.x = FALSE, sex.sp = TRUE)
#View(qtl::map2table(test.map))


### FldTrial.lm -- this fucntion is actually better, replace the other one with this
#Xlm <- function(y, X){
#  n <- nrow(X)
#  k <- ncol(X)
#  
#  inv.tX.X <- solve(crossprod(X))
#  #if(!exists("inv.tX.X")) inv.tX.X <- MASS::ginv(t(X) %*% X) ## If chol fails then use general inverse of MASS package
#  
#  beta <- matrix(inv.tX.X %*% t(X) %*% y, ncol = 1, dimnames = list(colnames(X), NULL))
#  
#  yhat <- X %*% beta
#  e <- as.vector(y - yhat)
#  
#  sigma.samp <- (t(e) %*% e) / (n-k)
#  varB <- sigma.samp %x% inv.tX.X
#  se.B <- sqrt(sigma.samp) %x% sqrt(diag(inv.tX.X))
#  rownames(varB) <- colnames(varB) <- rownames(beta)
#  
#  return(list(betahat=beta, resids=e, varbetahat= diag(varB), stderrbetahat=se.B, dfErr=(n-k)))
#}




## G.CV must have individuals as rows and markers as columns... so 100 individuals with 500 markers would be a 100x500 matrix
## y.CV must be a numeric vector representing a continuous variable


#'
#' @rdname internal
#' 
XValidate_nonInd <- function(y.CV=NULL, G.CV=NULL, models.CV=NULL, frac.train.CV=NULL, nCV.iter.CV=NULL, burnIn.CV=NULL, nIter.CV=NULL){
  
  gc(verbose = F) ## Close unused connections
  #con.path <- getwd() ## BGLR will write temp files to the wd
  
  non.BGLR <- models.CV[models.CV %in% c("rrBLUP")]
  BGLR <- models.CV[models.CV %in% c("BayesA", "BayesB", "BayesC", "BL", "BRR")]
  
  for(i in 1:nCV.iter.CV){
    if(i==1){
      rrBLUP.cv <- c()
      BayesA.cv <- c()
      BayesB.cv <- c()
      BayesC.cv <- c()
      BL.cv <- c()
      BRR.cv <- c()
    }
    
    TP.sample <- sample(x=1:length(y.CV), size=round(frac.train.CV*length(y.CV)), replace=F)
    VP.sample <- setdiff(1:length(y.CV), TP.sample)
    
    TP.G <- G.CV[TP.sample, ]
    TP.y <- y.CV[TP.sample]
    
    VP.G <- G.CV[VP.sample, ]
    VP.y <- y.CV[VP.sample]
    
    ### non-BGLR models below
    if("rrBLUP" %in% non.BGLR) {RR.pred <- rrBLUP::kinship.BLUP(y=TP.y, G.train=TP.G, G.pred=VP.G, K.method="RR"); rrBLUP.cv[i] <- cor(VP.y, RR.pred$g.pred, use="pairwise.complete.obs")}
    
    ### models from BGLR package
    ## Bayesian A
    if("BayesA" %in% BGLR){
      tryCatch({
        BayesA.fit <- BGLR::BGLR(y=TP.y, ETA=list(list(X=TP.G, model="BayesA")), verbose=F, nIter=1500, burnIn=500)
      }, warning=function(w){
        gc(verbose = F)
        BayesA.fit <- NULL
      }, error=function(e){
        gc(verbose = F)
        BayesA.fit <- NULL
      }
      ); gc(verbose = F)
      
      if(!is.null(BayesA.fit)){
        mkr.effs <- as.numeric(BayesA.fit$ETA[[1]]$b); BayesA.cv[i] <- cor((VP.G %*% mkr.effs), VP.y, use="pairwise.complete.obs")
      } else{BayesA.cv[i] <- NA}
    }
    
    ## Bayesian B 
    if("BayesB" %in% BGLR){
      tryCatch({
        BayesB.fit <- BGLR::BGLR(y=TP.y, ETA=list(list(X=TP.G, model="BayesB")), verbose=F, nIter=1500, burnIn=1000)
      }, warning=function(w){
        gc(verbose = F)
        BayesB.fit <- NULL
      }, error=function(e){
        gc(verbose = )
        BayesB.fit <- NULL
      }
      ); gc(verbose = F)
      
      if(!is.null(BayesB.fit)){
        mkr.effs <- as.numeric(BayesB.fit$ETA[[1]]$b)
        BayesB.cv[i] <- cor((VP.G %*% mkr.effs), VP.y, use="pairwise.complete.obs")
      } else{BayesB.cv[i] <- NA}
    }
    
    ## Bayesian C
    if("BayesC" %in% BGLR){
      tryCatch({
        BayesC.fit <- BGLR::BGLR(y=TP.y, ETA=list(list(X=TP.G, model="BayesC")), verbose=F, nIter=1500, burnIn=1000)
      }, warning=function(w){
        gc(verbose = F)
        BayesC.fit <- NULL
      }, error=function(e){
        gc(verbose = F)
        BayesC.fit <- NULL
      }
      ); gc(verbose = F)
      
      if(!is.null(BayesC.fit)){
        mkr.effs <- as.numeric(BayesC.fit$ETA[[1]]$b)
        BayesC.cv[i] <- cor((VP.G %*% mkr.effs), VP.y, use="pairwise.complete.obs")
      } else{BayesC.cv[i] <- NA}
    }
    
    ## Bayesian LASSO
    if("BL" %in% BGLR){
      tryCatch({
        BL.fit <- BGLR::BGLR(y=TP.y, ETA=list(list(X=TP.G, model="BL")), verbose=F, nIter=1500, burnIn=1000)
      }, warning=function(w){
        gc(verbose = F)
        BL.fit <- NULL
      }, error=function(e){
        gc(verbose = F)
        BL.fit <- NULL
      }
      ); gc(verbose = F)
      
      if(!is.null(BL.fit)){
        mkr.effs <- as.numeric(BL.fit$ETA[[1]]$b)
        BL.cv[i] <- cor((VP.G %*% mkr.effs), VP.y, use="pairwise.complete.obs")
      } else{BL.cv[i] <- NA}
    }
    
    ### Bayesian Ridge Reg.
    if("BRR" %in% BGLR){
      tryCatch({
        BRR.fit <- BGLR::BGLR(y=TP.y, ETA=list(list(X=TP.G, model="BRR")), verbose=F, nIter=1500, burnIn=1000)
      }, warning=function(w){
        gc(verbose = F)
        BRR.fit <- NULL
      }, error=function(e){
        gc(verbose = F)
        BRR.fit <- NULL
      }
      ); gc(verbose = F)
      
      if(!is.null(BRR.fit)){
        mkr.effs <- as.numeric(BRR.fit$ETA[[1]]$b)
        BRR.cv[i] <- cor((VP.G %*% mkr.effs), VP.y, use="pairwise.complete.obs")
      } else{BRR.cv[i] <- NA}
    }
    
  }; gc(verbose = F) ## End of ith iteration
  
  cv.list <- paste(models.CV, "cv", sep=".")
  
  for(l in 1:length(cv.list)){ ## Tried using sapply() across cv.list using the "get" function, but did not work within PopVar function
    toget <- cv.list[l]
    if(l == 1) cvs <- get(toget)
    if(l > 1) cvs <- cbind(cvs, get(toget))
  }
  
  if(length(models.CV) == 1){
    bad.models <- F
    if(length(which(is.na(cvs))) > 0.025*nCV.iter.CV) bad.models <- T
  }
  
  if(length(models.CV) > 1){
    bad.models <- apply(cvs, 2, function(X){if(length(which(is.na(X))) > 0.025*nCV.iter.CV){ ## if more than 2.5% of the iterations resulted in an error 
      return(T)
    }else{return(F)}
    })
  }
  
  if(length(models.CV) > 1) {CV.results <- data.frame(Model=models.CV, r_avg=apply(cvs, 2, mean, na.rm=T), r_sd=apply(cvs, 2, sd, na.rm=T)) ; rownames(CV.results) <- NULL}
  if(length(models.CV) == 1) {CV.results <- data.frame(Model=models.CV, r_avg=mean(cvs), r_sd=sd(cvs)) ; rownames(CV.results) <- NULL}
  
  if(length(which(bad.models)) > 0){
    if(length(models.CV) == 1 | (length(which(bad.models)) == length(models.CV))) stop("All model(s) tested was/were removed due to excessive negative values of nu being returned by BGLR::BGLR")
    CV.results <- CV.results[-which(bad.models), ]
    warning(paste("Model(s)", models.CV[which(bad.models)], "was/were removed due to excessive negative values of nu being returned by BGLR::BGLR."))
  }
  
  #CV.lists <- as.data.frame(t(rbind(as.character(models.CV), matrix(c(rrBLUP.cv, BayesA.cv, BayesB.cv, BayesC.cv, BL.cv, BRR.cv), ncol=length(models.CV)))))
  
  return(list(CV.summary = CV.results)) #, iter.CV = CV.lists))
  
} ## End of XValidate function



## G.CV must have individuals as rows and markers as columns... so 100 individuals with 500 markers would be a 100x500 matrix
## y.CV must be a numeric vector representing a continuous variable

#'
#' @rdname internal
#' 
XValidate_Ind <- function(y.CV=NULL, G.CV=NULL, models.CV=NULL, nFold.CV=NULL, nFold.CV.reps=NULL, burnIn.CV=NULL, nIter.CV=NULL){
  
  gc(verbose = F) ## Close unused connections
  #con.path <- getwd() ## BGLR will write temp files to the wd
  
  non.BGLR <- models.CV[models.CV %in% c("rrBLUP")]
  BGLR <- models.CV[models.CV %in% c("BayesA", "BayesB", "BayesC", "BL", "BRR")]
  if(nFold.CV >= length(y.CV)) stop("nFold too large given the TP size")
  
  for(j in 1:nFold.CV.reps){
    for(i in 1:nFold.CV){
      
      if(i==1){
        rrBLUP.cv <- c()
        BayesA.cv <- c()
        BayesB.cv <- c()
        BayesC.cv <- c()
        BL.cv <- c()
        BRR.cv <- c()
        
        for.VP.mat <- 1:length(y.CV)
        while(length(for.VP.mat) %% nFold.CV != 0) for.VP.mat <- c(for.VP.mat, NA)
        VP.mat <- matrix(for.VP.mat[order(sample(1:length(for.VP.mat), replace = F))], nrow = nFold.CV, byrow = T)
      }
      
      VP.sample <- as.numeric(VP.mat[i, ]); VP.sample <- VP.sample[which(!is.na(VP.sample))]
      TP.sample <- setdiff(1:length(y.CV), VP.sample)
      
      TP.G <- G.CV[TP.sample, ]
      TP.y <- y.CV[TP.sample]
      
      VP.G <- G.CV[VP.sample, ]
      VP.y <- y.CV[VP.sample]
      
      ### non-BGLR models below
      if("rrBLUP" %in% non.BGLR) {RR.pred <- rrBLUP::kinship.BLUP(y=TP.y, G.train=TP.G, G.pred=VP.G, K.method="RR"); rrBLUP.cv[i] <- cor(VP.y, RR.pred$g.pred, use="pairwise.complete.obs")}
      
      ### models from BGLR package
      ## Bayesian A
      if("BayesA" %in% BGLR){
        tryCatch({
          BayesA.fit <- BGLR::BGLR(y=TP.y, ETA=list(list(X=TP.G, model="BayesA")), verbose=F, nIter=1500, burnIn=1000)
        }, error=function(e){
          gc(verbose=F)
          BayesA.fit <- NULL
        }, warning=function(w){
          gc(verbose=F)
          BayesA.fit <- NULL
        }); gc(verbose=F)
        
        if(!is.null(BayesA.fit)){
          mkr.effs <- as.numeric(BayesA.fit$ETA[[1]]$b); BayesA.cv[i] <- cor((VP.G %*% mkr.effs), VP.y, use="pairwise.complete.obs")
        } else{BayesA.cv[i] <- NA}
      }
      
      
      ## Bayesian B 
      if("BayesB" %in% BGLR){
        tryCatch({
          BayesB.fit <- BGLR::BGLR(y=TP.y, ETA=list(list(X=TP.G, model="BayesB")), verbose=F, nIter=1500, burnIn=1000)
        }, error=function(e){
          gc(verbose=F)
          BayesB.fit <- NULL
        }, warning=function(w){
          gc(verbose=F)
          BayesB.fit <- NULL
        }); gc(verbose=F)
        
        if(!is.null(BayesB.fit)){
          mkr.effs <- as.numeric(BayesB.fit$ETA[[1]]$b)
          BayesB.cv[i] <- cor((VP.G %*% mkr.effs), VP.y, use="pairwise.complete.obs")
        } else{BayesB.cv[i] <- NA}
      }
      
      
      ## Bayesian C
      if("BayesC" %in% BGLR){
        tryCatch({
          BayesC.fit <- BGLR::BGLR(y=TP.y, ETA=list(list(X=TP.G, model="BayesC")), verbose=F, nIter=1500, burnIn=1000)
        }, error=function(e){
          gc(verbose=F)
          BayesC.fit <- NULL
        }, warning=function(w){
          gc(verbose=F)
          BayesC.fit <- NULL
        }); gc(verbose=F)
        
        if(!is.null(BayesC.fit)){
          mkr.effs <- as.numeric(BayesC.fit$ETA[[1]]$b)
          BayesC.cv[i] <- cor((VP.G %*% mkr.effs), VP.y, use="pairwise.complete.obs")
        } else{BayesC.cv[i] <- NA}
      }
      
      
      ## Bayesian LASSO
      if("BL" %in% BGLR){
        tryCatch({
          BL.fit <- BGLR::BGLR(y=TP.y, ETA=list(list(X=TP.G, model="BL")), verbose=F, nIter=1500, burnIn=1000)
        }, error=function(e){
          gc(verbose=F)
          BL.fit <- NULL
        }, warning=function(w){
          gc(verbose=F)
          BL.fit <- NULL
        }); gc(verbose=F)
        
        if(!is.null(BL.fit)){
          mkr.effs <- as.numeric(BL.fit$ETA[[1]]$b)
          BL.cv[i] <- cor((VP.G %*% mkr.effs), VP.y, use="pairwise.complete.obs")
        } else{BL.cv[i] <- NA}        
      }
      
      
      ### Bayesian Ridge Reg.
      if("BRR" %in% BGLR){
        tryCatch({
          BRR.fit <- BGLR::BGLR(y=TP.y, ETA=list(list(X=TP.G, model="BRR")), verbose=F, nIter=1500, burnIn=1000)
        }, error=function(e){
          gc(verbose=F)
          BRR.fit <- NULL
        }, warning=function(w){
          gc(verbose=F)
          BRR.fit <- NULL
        }); gc(verbose=F)
        
        if(!is.null(BRR.fit)){
          mkr.effs <- as.numeric(BRR.fit$ETA[[1]]$b)
          BRR.cv[i] <- cor((VP.G %*% mkr.effs), VP.y, use="pairwise.complete.obs")
        } else{BRR.cv[i] <- NA}
      }
      
    }; gc(verbose=F) ## End of i-fold iteration
    
    if(j==1) cv.list <- paste(models.CV, "cv", sep=".")
    
    for(l in 1:length(cv.list)){ ## Tried using sapply() across cv.list using the "get" function, but did not work within PopVar function
      toget <- cv.list[l]
      if(l == 1) cvs <- get(toget)
      if(l > 1) cvs <- cbind(cvs, get(toget))
    }
    
    if(j == 1) cvs.all <- cvs
    if(j > 1) cvs.all <- rbind(cvs.all, cvs)
    
  }
  
  if(length(models.CV) == 1){
    bad.models <- F
    if(length(which(is.na(cvs.all))) > 0.025*(nFold.CV*nFold.CV.reps)) bad.models <- T
  }
  
  if(length(models.CV) > 1){
    bad.models <- apply(cvs.all, 2, function(X){if(length(which(is.na(X))) > 0.025*(nFold.CV*nFold.CV.reps)){ ## if more than 2.5% of the iterations resulted in an error 
      return(T)
    }else{return(F)}
    })
  }
  
  if(length(models.CV) > 1) {CV.results <- data.frame(Model=models.CV, r_avg=apply(cvs.all, 2, mean, na.rm=T), r_sd=apply(cvs.all, 2, sd, na.rm=T)) ; rownames(CV.results) <- NULL}
  if(length(models.CV) == 1) {CV.results <- data.frame(Model=models.CV, r_avg=mean(cvs.all), r_sd=sd(cvs.all)) ; rownames(CV.results) <- NULL}
  
  if(length(which(bad.models)) > 0){
    if(length(models.CV) == 1 | (length(which(bad.models)) == length(models.CV))) stop("All model(s) tested was/were removed due to excessive negative values of nu being returned by BGLR::BGLR")
    CV.results <- CV.results[-which(bad.models), ]
    warning(paste("Model(s)", models.CV[which(bad.models)], "was/were removed due to excessive negative values of nu being returned by BGLR::BGLR."))
  }
  
  #CV.lists <- as.data.frame(t(rbind(as.character(models.CV), matrix(c(rrBLUP.cv, BayesA.cv, BayesB.cv, BayesC.cv, BL.cv, BRR.cv), ncol=length(models.CV)))))
  return(list(CV.summary = CV.results))#, iter.CV = CV.lists))
  
} ## End of XValidate function





# Function to calculate marker effects
# 
# Allows other arguments to be passed
# 

#'
#' @rdname internal
#' 
calc_marker_effects <- function(M, y.df, models = c("rrBLUP", "BayesA", "BayesB","BayesC", "BL", "BRR"), nIter, burnIn) {
  
  models <- match.arg(models)
  
  # Error out if model != rrBLUP and nIter & burnIn are missing
  if (models != "rrBLUP") {
    if (missing(nIter) | missing(burnIn)) stop ("You must provide the arguments 'nIter' and 'burnIn' for that model choice.")
    
  }
  
  # Determine the function to use for marker effect estimation
  if (models == "rrBLUP") {
    cme <- function(M, y) {
      fit <- mixed.solve(y = y, Z = M, method = "REML")
      # Return marker effects and the grand mean
      list(effects = as.matrix(fit$u), grand_mean = fit$beta)
    }
  } else {
    cme <- function(M, y) {
      suppressWarnings(fit <- BGLR(y = y, ETA = list(list(X = M, model = models)), verbose = FALSE, nIter = nIter, burnIn = burnIn))
      list(effects = as.matrix(fit$ETA[[1]]$b), grand_mean = fit$mu)
    }
  }
    
  ## Calculate marker effects for each trait
  me_out <- lapply(X = y.df, FUN = cme, M = M)
  
  # Clean up files if models != rrBLUP
  if (models != "rrBLUP") {
    invisible(file.remove(list.files(path = ".", pattern = paste0("mu|varE|", models), full.names = TRUE)))
  }
 
  # Return me_out
  return(me_out)
   
}


#' Filter SNP marker genotypes
#'
#' @param geno A data.frame of genotype calls, where the first three columns are "marker," "chrom,", and "position."
#' Subsequent columns are genotype calls for individuals.
#' @param max.ind.miss The maximum rate of missingness to retain an individual
#' @param max.mar.miss The maximum rate of missingness to retain an marker
#' @param min.mac The minimum minor allele count to retain a marker.
#' @param max.r2 The maximum squared correlation between any pair of markers. The \code{r^2} is used directly in filtering. For each pair of markers exceeding this threshold, the one with the greater minor allele frequency is retained. Set \code{max.r2 = 1} for no filtering.
#'
#' @return
#' Variable of class \code{\link{geno}}.
#'
#' @examples
#' # Use data from the PopVar package
#' data("cranberry_geno", package = "PopVar")
#'
#' # Remove cM, rename
#' geno <- cranberry_geno[-3]
#' names(geno)[1:3] <- c("marker", "chrom", "position")
#'
#' # Filter
#' geno_in <- filter_geno(geno = geno, max.ind.miss = 0.2, max.mar.miss = 0.2, min.mac = 5, max.r2 = 1)
#'
filter_geno <- function(geno, max.ind.miss = 0.5, max.mar.miss = 0.5, min.mac = 5, max.r2 = 1) {
  
  stopifnot(is.data.frame(geno))
  cols.match <- names(geno)[1:3] == c("marker", "chrom", "position")
  if (any(!cols.match)) {
    stop("The first 3 columns of 'geno' should be 'marker', 'chrom', and 'position'.")
  }
  
  # Error checking
  stopifnot(max.mar.miss >= 0 & max.mar.miss <= 1)
  stopifnot(max.ind.miss >= 0 & max.ind.miss <= 1)
  stopifnot(min.mac >= 0)
  stopifnot(max.r2 >= 0 & max.r2 <= 1)
  
  geno_orig <- geno
  
  geno <- as.matrix(geno_orig[,-1:-3])
  row.names(geno) <- all.marks <- geno_orig$marker
  all.marks <- row.names(geno)
  n <- ncol(geno)
  p <- nrow(geno)
  ix <- which(apply(geno, 1, sd, na.rm = T) > 0)
  nd <- nrow(geno) - length(ix)
  cat(sub("X", nd, "Removed X markers without genetic variance\n"))
  geno <- geno[ix, , drop = FALSE]
  
  mar_miss <- rowMeans(is.na(geno))
  iu <- which(mar_miss <= max.mar.miss)
  nd2 <- nrow(geno) - length(iu)
  cat(sub("X", nd2, "Removed X markers due to missing data\n"))
  geno <- geno[iu, , drop = FALSE]
  
  ind_miss <- colMeans(is.na(geno))
  iv <- which(ind_miss <= max.ind.miss)
  ni <- ncol(geno) - length(iv)
  cat(sub("X", ni, "Removed X individuals due to missing data\n"))
  geno <- geno[, iv, drop = FALSE]
  
  # Calculate maf
  maf <- calc_maf(x = t(geno) - 1, check.matrix = FALSE)
  # How many markers
  iw <- which(maf >= min.mac / ncol(geno))
  nm <- nrow(geno) - length(iw)
  cat(sub("X", nm, "Removed X markers due to low minor allele count\n"))
  geno <- geno[iw, , drop = FALSE]
  
  # Compute pairwise r2 - only if max.r2 < 1
  if (max.r2 < 1) {
    # Impute genotypes with the mean, compute the correlation
    geno1 <- apply(X = geno, MARGIN = 1, FUN = function(snp) {
      snp[is.na(snp)] <- mean(snp, na.rm = TRUE)
      snp
    })
    geno1 <- geno1 - 1
    
    mar_r2 <- cor(x = geno1)^2
    
    # Run SNP pruning
    pruned_ans <- prune_LD(x = geno1, cor.mat = mar_r2, r2.max = max.r2, check.matrix = FALSE)
    # Report
    mars_kept <- colnames(pruned_ans)
    nr <- nrow(geno) - length(mars_kept)
    cat(sub("X", nr, "Removed X markers due to high pairwise LD\n"))
    
    geno <- geno[mars_kept, , drop = FALSE]
    
  }
  
  # Get the marker information back
  snp_info <- geno_orig[c(1:3)]
  snp_info <- snp_info[snp_info$marker %in% row.names(geno), ]
  geno_out <- cbind(snp_info[match(x = row.names(geno), table = snp_info$marker), ], geno)
  row.names(geno_out) <- NULL
  return(geno_out)
  
}


#' Prune SNPs based on linkage disequilibrium
#'
#' @param x An n x p marker matrix of n individuals and p markers coded as -1, 0, and 1
#' @param cor.mat A similarity matrix between SNPs. No calculation of LD will be made if this matrix is passed.
#' for homozygous alternate, heterozygous, and homozygous reference.
#' @param r2max The maximum acceptable value of r2 (i.e. LD) or similarity between any two markers.
#' @param check.matrix Logical. Should the marker matrix 'x' be checked? Use check.matrix = FALSE
#' if 'x' contains imputed decimal genotypes.
#'
#' @return
#' A marker matrix without markers in high LD
#' 
prune_LD <- function(x, cor.mat, r2.max = 0.80, check.matrix = TRUE) {
  
  ## Error checking
  # Check the marker matrix
  stopifnot(is.logical(check.matrix))
  if (check.matrix) stopifnot(check_marker_matrix(x))
  if (!missing(cor.mat)) {
    stopifnot(is.matrix(cor.mat))
  }
  # check r2max
  stopifnot(r2.max >= 0 & r2.max <= 1)
  
  # calculate minor allele frequency
  maf <- calc_maf(x, check.matrix = check.matrix)
  
  if (missing(cor.mat)) {
    
    # calculate the correlation between all markers; square it
    if (any(is.na(x))) {
      all_marker_r <- cor(x, use = "pairwise.complete.obs")^2
      
    } else {
      all_marker_r <- cor(x)^2
    }
    
  } else {
    all_marker_r <- cor.mat
    
  }
  
  
  
  # Set the lower half (including diagonal) to NA
  all_marker_r1 <- all_marker_r
  all_marker_r1[lower.tri(all_marker_r1, diag = TRUE)] <- NA
  
  # Get a matrix of those entries that are elevated
  elevated_marker_r <- all_marker_r1 > r2.max
  # Get the coordinates of those entries
  which_elevated <- which(x = elevated_marker_r, arr.ind = TRUE)
  
  markers_remove <- character()
  # While loop
  i = 1
  while(nrow(which_elevated) > 0) {
    
    # Subset the first row
    coords <- which_elevated[1,]
    # Extract the coordinate
    r2_coord <- all_marker_r1[coords[1], coords[2], drop = FALSE]
    
    # marker pair
    markers <- unlist(dimnames(r2_coord))
    # Identify the marker with higher MAF
    higher_maf <- which.max(maf[markers])
    
    # Toss that marker
    marker_remove <- names(higher_maf)
    markers_remove[i] <- marker_remove
    
    # Find the row/col containing this marker
    row_remove <- col_remove <- which(row.names(all_marker_r1) == marker_remove)
    
    which_elevated <- subset.matrix(which_elevated, which_elevated[,"row"] != row_remove & which_elevated[,"col"] != col_remove)
    
    # advance i
    i <- i + 1
    
  }
  
  # Remove the markers from the marker matrix
  cols_keep <- setdiff(seq_len(ncol(x)), which(colnames(x) %in% markers_remove))
  x[,cols_keep,drop = FALSE]
  
}


#' Calculate minor allele frequency
#'
#' @param x An n x p marker matrix of n individuals and p markers coded as -1, 0, and 1
#' for homozygous alternate, heterozygous, and homozygous reference.
#' @param check.matrix Logical. Should the marker matrix 'x' be checked? Use check.matrix = FALSE
#' if 'x' contains imputed decimal genotypes.
#'
calc_maf <- function(x, check.matrix = TRUE) {
  
  ## Error checking
  stopifnot(is.logical(check.matrix))
  if (check.matrix) stopifnot(check_marker_matrix(x))
  
  # Filter for MAF
  af <- colMeans(x + 1, na.rm = TRUE) / 2
  maf <- pmin(af, 1 - af)
  
  return(maf)
  
}











