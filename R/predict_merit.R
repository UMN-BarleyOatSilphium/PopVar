#' Predict breeding values or genotypic values
#' 
#' 
#' @param data A \code{PopVar.data} object from \code{read_PopVar}.
#' @param models The list of models to run. Currently supports ...
#' @param method The cross-validation method. Either fractional (\code{"frac"}), k-fold (\code{"kfold"}), or leave-one-out \code{"loo"}.
#' @param p.test For fractional cross-validation, the proportion of individuals to designate as the test set per replicate. \code{1 - p.test} is the proportion of individuals used for training.
#' @param k For k-fold cross-validation, the number of folds.
#' @param reps The number of complete replicates of fractional or k-fold cross-validation.
#' @param group The name of a column in the phenotype data for grouped cross-validation (e.g. environment, location, etc.). For example to run leave-one-environment-out cross-validation, set \code{method = "loo"} and \code{group = "env"}.
#' 
#' 
#' @export
#' 
predict_merit <- function(data) {}
  