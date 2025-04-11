#' Predict marker genetic position using linear interpolation
#'
#' @description
#' Uses linear interpolation and local recombination rate information to predict
#' the unknown genetic position of markers based on their known physical positions.
#'
#' @param geno A data frame of marker allele dosages. See \italic{Details} for more information.
#' @param cMMb A data frame of the recombination rate information. See \italic{Details}.
#' 
#' @details
#' The first columns of the geno object are marker, chrom, and bp. Subsequent columns 
#' contain the allele dosage for individuals/clones, coded 0,1,2. While fractional values are allowed for genomewide prediction, 
#' these values are rounded when predicting the genetic variance within crosses.
#' 
#' The cMMb object must contains six columns: chrom, left_bp (the physical position of the left
#' flank of an interval), right_bp (the physical position of the right flank of an interval),
#' left_cM (the genetic position of the left flank of an interval), right_cM 
#' (the genetic position of the right flank of an interval), and cM_Mb (the recombination rate 
#' (in cM / Mb) within that interval).
#'
#' @export
#'
interpolate_genetic_position <- function(geno, cMMb) {
  # Error checking
  stopifnot(is.data.frame(geno))
  stopifnot(is.data.frame(cMMb))
  
  # Extract the map information
  cols_required <- c("marker", "chrom", "bp")
  if (any(!cols_required %in% colnames(geno))) {
    stop("The first three columns of 'filename' or 'geno' must be 'marker', 'chrom', and 'bp'.")
  }
  if (colnames(geno)[4] == "bp") {
    map <- geno[,1:4]
    geno <- as.matrix(geno[,-(1:4)])
  } else {
    map <- geno[,1:3]
    geno <- as.matrix(geno[,-(1:3)])
  }
  
  # Get the recombination rate information
  cols_required <- c("chrom", "left_bp", "right_bp", "left_cM", "right_cM", "cM_Mb")
  if (any(!cols_required %in% colnames(cMMb))) {
    stop("The first three columns of 'filename' or 'geno' must be 'marker', 'chrom', and 'bp'.")
  }
  cMMb1 <- cMMb[,cols_required]
  
  # All chromosomes in map in cmmb?
  chroms_map <- unique(map$chrom)
  chroms_cmmb <- unique(cMMb1$chrom)
  
  if (!all(chroms_map %in% chroms_cmmb)) stop ("The chromosomes in 'geno' are not all in 'cMMb'.")
  
  # Make sure the cMMb starts at zero
  cMMb_split <- split(cMMb1, cMMb1$chrom)
  cMMb_split1 <- cMMb_split
  for (i in seq_along(cMMb_split)) {
    cMMb_i <- cMMb_split[[i]]
    if (cMMb_i$left_bp[1] > 0) {
      head <- data.frame(chrom = unique(cMMb_i$chrom), left_bp = 0, right_bp = cMMb_i$left_bp[1], left_cM = 0, right_cM = cMMb_i$left_cM[1], cM_Mb = as.numeric(NA))
      head$cM_Mb <- (head$right_cM - head$left_cM) / ((head$right_bp - head$left_bp) / 1e6)
      cMMb_i1 <- rbind(head, cMMb1)
      
    }
    row.names(cMMb_i1) <- NULL
    cMMb_split1[[i]] <- cMMb_i1
    
  }
  
  cMMb2 <- do.call(rbind, cMMb_split1)
  row.names(cMMb2) <- NULL
  
  
  # Add a new column to map
  map1 <- map
  map1$cM <- as.numeric(NA)
  
  map1_split <- split(map1, map1$chrom)
  
  # Iterate over chromosomes
  for (i in seq_along(map1_split)) {
    map_i <- map1_split[[i]]
    cmmb_i <- cMMb2[cMMb2$chrom == unique(map_i$chrom), ]
    
    # Create a set of x-y coordinates
    n <- nrow(cmmb_i)
    cmmb_i_xy <- cmmb_i[,c("chrom", "left_bp", "left_cM")]
    colnames(cmmb_i_xy) <- c("chrom", "bp", "cM")
    cmmb_n <- cmmb_i[n,c("chrom", "right_bp", "right_cM")]
    colnames(cmmb_n) <- c("chrom", "bp", "cM")
    cmmb_i_xy <- rbind(cmmb_i_xy, cmmb_n)
    
    map_max_bp <- max(map_i$bp)
    if (map_max_bp > max(cmmb_i_xy$bp)) {
      cmmb_n[,"bp"] <- map_max_bp
      cmmb_i_xy <- rbind(cmmb_i_xy, cmmb_n)
    }
    
    # Interpolate
    interp_out <- approx(x = cmmb_i_xy$bp, y = cmmb_i_xy$cM, xout = map_i$bp, ties = "ordered")
    # Add the cM to the map
    map_i$cM <- interp_out$y
    
    map1_split[[i]] <- map_i
    
  }
  
  map2 <- do.call(rbind, map1_split)
  row.names(map2) <- NULL
  
  # Merge with the genotypes
  geno1 <- cbind(map2, geno)
  row.names(geno1) <- NULL
  
  
  # Return the adjusted genotype 
  return(geno1)
  
}
