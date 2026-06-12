#' @title The function of METEOR method
#' @description  METEOR is able to model multiple outcomes simultaneously with an initial set of candidate SNP instruments that are in high LD with each other and perform automated instrument selection to identify suitable SNPs to serve as instrumental variables.
#' @param Zscore the Zscore with k+1 rows and p columns of the SNP effect size matrix for the exposure and k outcomes
#' @param Sigma the LD matrix for the SNPs selected from the exposure can be obtained by using the weighted average LD matrix from k outcomes. If individual data is unavailable, the LD matrix can also be derived from a reference panel.
#' @param N the sample size of exposure and k outcomes GWASs
#' @param Omega_est the correlation matrix derived by LDSC
#' @param Omega_se the standard errors of elements of Omega derived by LDSC
#' @param Gibbsnumber the number of Gibbs sampling iterations with the default to be 1000
#' @param burninproportion  the proportion to burn in from Gibbs sampling iterations, with default to be 0.2  
#' @param lambda the tuning parameter used to ensure that the correlation matrix is invertible
#' @param pi_beta_shape the prior shape paramter for pi_beta with the default to be 0.5
#' @param pi_beta_scale the prior scale paramter for pi_beta with the default to be 4.5
#' @param pi_1_shape the prior shape paramter for pi_1 with the default to be 0.5
#' @param pi_1_scale the prior scale paramter for pi_1 with the default to be 1.5
#' @param pi_0_shape the prior shape paramter for pi_0 with the default to be 0.05
#' @param pi_0_scale the prior scale paramter for pi_0 with the default to be 9.95
#' @return A list of estimated parameters for the causal effect test 
#' \item{causal_effect}{The estimate of k causal effects}
#' \item{causal_pvalue_single}{The p values for the causal effects in single tests}
#' \item{causal_pvalue_overall}{The p value for the causal effects in overall test}
#' \item{causal_sd}{The standard deviation for the causal effects}
#' \item{causal_cov}{The covariance for the causal effects}
#' \item{sigmabeta}{The variance estimate for the SNP effect sizes on the exposure}
#' \item{sigmaeta}{The variance estimates for the horizontal pleiotropy effects}
#' \item{Omega}{The correlation matrix for one exposure and k outcomes}
#' \item{pi_beta}{The proportion of selected SNPs, which show non-zero effects on exposure}
#' \item{pi_1}{The proportion of selected SNPs showing horizontal pleiotropy}
#' \item{pi_0}{The proportion of non-selected SNPs showing horizontal pleiotropy}

METEOR<-function(Zscore,Sigma,N,Omega_est,Omega_se,Gibbsnumber=1000,burninproportion=0.2,lambda=0.1,pi_beta_shape=0.5,
pi_beta_scale=4.5,pi_1_shape=0.5,pi_1_scale=1.5,pi_0_shape=0.05,pi_0_scale=9.95){
  Sigma <- as.matrix(Sigma)
  Omega_est <- as.matrix(Omega_est)
  Omega_se <- as.matrix(Omega_se)
  pvalue<-2*(1-pnorm(abs(Omega_est/Omega_se)))
  k <- nrow(Omega_est)
  for (m1 in 1:(k-1)){
    for (m2 in (m1+1):k){
      if (pvalue[m1,m2]>0.05){
        Omega_est[m1,m2] <- 0
        Omega_est[m2,m1] <- 0 
      }
    }
  }
  Omega <- Omega_est

  Omega <- as.matrix(Omega)
  p = ncol(Sigma)
  k = nrow(Zscore)-1
  if (dim(Zscore)[2] == p){
    if (det(Omega)<=0){
      I = diag(nrow(Omega))
      Omega = Omega+lambda*I
    }
    re=METEOR_CPP(Zscore_in=Zscore,Sigma_in=Sigma,N_in=N,p_in = p,k_in=k,Omega_in=Omega,Gibbsnumberin=Gibbsnumber,
                 burninproportion=burninproportion,pi_beta_shape_in=pi_beta_shape,pi_beta_scale_in=pi_beta_scale,
                 pi_1_shape_in=pi_1_shape,pi_1_scale_in=pi_1_scale,pi_0_shape_in=pi_0_shape,pi_0_scale_in=pi_0_scale)
    qq1 <- re$alpha*(1/(re$alpha_sd*re$alpha_sd))*re$alpha
    pvalue_single <- t(pchisq(t(qq1), 1, lower.tail = F))
    if (k>1){
      Aeigen=eigen(re$alpha_cov)
      WW=diag(1/Aeigen$values)
      U=Aeigen$vectors
      inv_alpha_cov =U%*%WW%*%t(U)
      qq2 = re$alpha%*%inv_alpha_cov%*%t(re$alpha)
      pvalue_overall = pchisq(qq2, k, lower.tail = F)
    } else {
      pvalue_overall = pvalue_single
    }
    
    result=list()
    result$causal_effect=re$alpha
    result$causal_pvalue_single=pvalue_single
    result$causal_pvalue_overall=pvalue_overall
    result$causal_sd=re$alpha_sd
    result$causal_cov=re$alpha_cov
    result$sigmabeta=re$sigma2beta
    result$sigmaeta=re$sigma2eta
    result$Omega=Omega
    result$pi_beta=re$pi_beta
    result$pi_1=re$pi_1
    result$pi_0=re$pi_0
    return(result)
  } else {
    print("ERROR: The dimensions of Zscore and Sigma are not matched. Zscore must be a matrix with k+1 rows and p columns, Sigma must be a matrix with p rows and p columns, where 'p' is the number of SNPs and 'k' is the number of outcomes")
  }
}
  