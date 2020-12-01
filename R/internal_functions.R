## Group of functions that PopVar uses internally



tails <- function(GEBVs, tail.p){ #Calculates means of tails; set tail.p to the proportion of the populaiton you want to take the mean of, default is 10%
  u.top <- mean(GEBVs[which(GEBVs >= quantile(GEBVs, 1-tail.p))], na.rm=T)
  u.bot <- mean(GEBVs[which(GEBVs <= quantile(GEBVs, tail.p))], na.rm=T)
  
  return(rbind(u.top, u.bot))
}

maf.filt <- function(G){
  G_noNA <- sum(!is.na(G), na.rm = TRUE)
  freq1 <- sum(G == 1, na.rm = TRUE) / G_noNA #+ .5*(sum(G == 0, na.rm = TRUE) / G_noNA)
  min(freq1, 1 - freq1)
}

cor.filt <- function(G, cor.cutoff){
  G4cor <- t(G)
  G4cor[which(is.na(G4cor))] <- 0
  G.cor <- cor(G4cor)
  G.cor[lower.tri(G.cor, diag = TRUE)] <- NA
  dup.mat <- which(G.cor > cor.cutoff, arr.ind = TRUE)
  dup.list <- as.numeric(dup.mat)
  two.dups <- as.numeric(names(which(table(dup.list) > 1)))
  dup.mat <- dup.mat[-sapply(two.dups, function(i) which(dup.mat == i, arr.ind = TRUE)[-1,1]), ]
  dup.mat[,2]
}

#test.map <- qtl::sim.map(len = c(50,50), n.mar = c(100,100), anchor.tel = FALSE, include.x = FALSE, sex.sp = TRUE)
#View(qtl::map2table(test.map))


mkr.sel <- function(){
  G.sig <- t(G_GWAS[,-c(1:3)]); colnames(G.sig) <- G_GWAS[,1]
  sig.mat <- data.frame(y_GWAS[,2], G.sig[,G.markers %in% sig.mkrs]); colnames(sig.mat)[1] <- t
  
  lm.form <- as.formula(paste(t, "~", paste(sig.mkrs, collapse = "+")))
  step.out <- step(lm(lm.form, data = sig.mat), scope = as.formula(paste("~", paste(sig.mkrs, collapse = "+"))), direction = "both", k = log(nrow(y_GWAS)), verbose = FALSE)
  list(final.mkrs = names(step.out$coefficients)[names(step.out$coefficients) %in% sig.mkrs], mkr.LD <- cor(sig.mat[,-1]))
}



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

XValidate.nonInd <- function(y.CV=NULL, G.CV=NULL, models.CV=NULL, frac.train.CV=NULL, nCV.iter.CV=NULL, burnIn.CV=NULL, nIter.CV=NULL){
  
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

XValidate.Ind <- function(y.CV=NULL, G.CV=NULL, models.CV=NULL, nFold.CV=NULL, nFold.CV.reps=NULL, burnIn.CV=NULL, nIter.CV=NULL){
  
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



## G.CV must be an n x m (marker) matrix
## y.CV must be a numeric vector representing a continuous variable


#y.CV = y_TP
#X.CV = X_TP
#G.CV = G_TP
#nFold.reps= 5
#model = "RR"
#CV.iter = 5000
#CV.burn = 1000
##


Xval <- function(model, y.CV, X.CV, G.CV){
  
  if(model == "BA"){iv = TRUE; pi = 0; de=FALSE}
  if(model == "BB"){iv = TRUE; pi = .95; de=FALSE}
  if(model == "BC"){iv = FALSE; pi = .95; de=FALSE}
  if(model == "BL"){iv = TRUE; pi = 0; de=TRUE}
  accL <- list()
  
  if(length(models) > 1){
    for(n in 1:nFold.reps){
      acc <- c()
      for(i in 1:nFolds){
        y_train <- y.CV[sets.train[[n]][[i]]]
        y_pred <- y.CV[sets.pred[[n]][[i]]]
        
        X_train <- X.CV[sets.train[[n]][[i]],]
        X_pred <- X.CV[sets.pred[[n]][[i]],]
        
        G_train <- G.CV[sets.train[[n]][[i]],]
        G_pred <- G.CV[sets.pred[[n]][[i]],]
        
        if(model %in% c("RR", "dnRR")){
          RRB.out <- rrBLUP::mixed.solve(y_train, X = X_train, Z = G_train)
          b <- RRB.out$beta; u <- RRB.out$u
          acc <- c(acc, cor(G_pred %*% u + as.matrix(X_pred) %*% b, y_pred))
        } else{
          bWGR.out <- bWGR::wgr(y = y_train, X = G_train, iv = iv, pi = pi, de = de, it = CV.iter, bi = CV.burn, verb = verbose)
          u <- bWGR.out$b
          acc <- c(acc, cor(G_pred %*% u , y_pred))
        }
      }
      accL[[n]] <- acc
    }  
  } else accL <- "CV not performed per user choice"
  
  ## Estimate mkr effs using full TP
  if(model %in% c("RR", "dnRR")){
    RRB.out <- rrBLUP::mixed.solve(y.CV, X = X.CV, Z = G.CV)
    b <- RRB.out$beta; u <- RRB.out$u
  } else{
    bWGR.out <- bWGR::wgr(y = y_train, X = G_train, iv = iv, pi = pi, de = de, it = CV.iter, bi = CV.burn, verb = verbose)
    u <- bWGR.out$b; b <- bWGR.out$mu
  }
  
  return(list(r = accL, u = u, b = b))
}



# Function to calculate marker effects
# 
# Allows other arguments to be passed
# 
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





