#' @title Estimate of correlation matrix without parallel.
#' @description Estimate of correlation matrix.
#'
#' @param data_list  a list of GWAS-summary-level data for exposure and k outcomes
#' @param ldscore.dir specify the path to the LD score files.
#' 
#' @export
#' @return A list of estimated parameters for the correlation matrix
#' \item{Omega}{the estimate of Omega }
#' \item{Omega.se}{the estimated matrix consists of the standard errors}
#' @importFrom parallel detectCores

Omega_est_nopar = function(data_list=summarydata,
                     ldscore.dir = ldscore.dir){
  trait_num = length(data_list)
  comb <- combn(seq_along(data_list), 2, simplify = FALSE)
  covariances <- list()
  for (i in 1:length(comb)){
    idx <- comb[[i]]
    data1 <- data_list[[idx[1]]]
    data2 <- data_list[[idx[2]]]
    paras = est_SS(dat1 = data1,
                   dat2 = data2,
                   trait1.name = "exp",
                   trait2.name = "out",
                   ldscore.dir = ldscore.dir)
    covariances[[i]] <- cbind(paras$C,paras$C.se)
  }
    
  names(covariances) <- sapply(comb, function(idx) paste0("Cov_", idx[1], "_", idx[2]))
  
  cov_matrix <- matrix(0, nrow = trait_num, ncol = trait_num)
  cov_matrix_se <- matrix(0,nrow = trait_num, ncol = trait_num)
  
  for (i in 1:trait_num){
    aa <- 0
    bb <- 0
    cc <- 0
    for (j in seq_along(comb)){
      idx <- comb[[j]]
      if (idx[1]==i){
        aa <- aa+covariances[[j]][1,1]
        cc <- cc+covariances[[j]][1,3]
        bb <- bb+1
      }
      if (idx[2]==i){
        aa <- aa+covariances[[j]][2,2]
        cc <- cc+covariances[[j]][2,4]
        bb <- bb+1
      }
    }
    cov_matrix[i,i] <- aa/bb
    cov_matrix_se[i,i] <- cc/bb
  }
  
  for (i in seq_along(comb)) {
    idx <- comb[[i]]
    cov_matrix[idx[1], idx[2]] <- covariances[[i]][1,2]
    cov_matrix[idx[2], idx[1]] <- covariances[[i]][1,2]
    cov_matrix_se[idx[1], idx[2]] <- covariances[[i]][1,4]
    cov_matrix_se[idx[2], idx[1]] <- covariances[[i]][1,4]
  }
  result_cov <- list()
  result_cov$Omega <- cov_matrix
  result_cov$Omega_se <- cov_matrix_se
  return(result_cov)
}

