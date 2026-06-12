---
layout: page
title: Example Analysis
description: ~
---
We consider a simple situation with two outcomes as an example for METEOR. Before running the tutorial, make sure that the MAPLE package is successfully installed. Please see the [link](https://Liye222.github.io/METEOR/documentation/02_installation.html) for the installation instructions. The example data for the tutorial can be downloaded in this [link](https://github.com/Liye222/METEOR/tree/main/example).

## A simulated data example

We conducted an alternative simulation in which exposure  exerts a non-zero effect on the first outcome while having no effect on the second outcome. This simulation was based on the realistic genotypes from UK Biobank, with sample sizes of n<sub>1</sub> = n<sub>21</sub> =  n<sub>22</sub> = 50,000 on chr 1. 

### Step 1: Estimation of correlation matrix

The function `Omega_est` or `Omega_est_nopar` can estimate the correlation matrix used to account for sample structure (e.g., population stratification, sample overlap and any correlations among them).

```r
library(Rcpp)
library(RcppArmadillo)
library(RcppDist)
library(magrittr)
library(data.table)
library(METEOR)
#load the summary data for exposure and two outcomes
load(file=paste0("./exposure.rda")) ##file name: exp
dat_x <- exp[,c("SNP","b","se","frq_A1","A1","A2","P","N")] 
load(file=paste0("./outcome1.rda")) ##file name: out1
dat_y1 <- out1[,c("SNP","b","se","frq_A1","A1","A2","P","N")] 
load(file=paste0("./outcome2.rda")) ##file name: out2
dat_y2 <- out2[,c("SNP","b","se","frq_A1","A1","A2","P","N")]

summarydata <- list(dat_x,dat_y1,dat_y2)
```
Users can use two functions, `Omega_est` or `Omega_est_nopar`, to estimate $\pmb\Omega$. The former utilizes `parallel`, which accelerates the computing time but requires more cores. The latter does not deed multiple cores, but its computing speed is relatively slow.

Function 1: `Omega_est`

```r
library(parallel)
Omega <- Omega_est(data_list=summarydata,
                   ldscore.dir = "./eur_w_ld_chr",
                   nCores=NA,
                   system_used="linux")
```

Function 2: `Omega_est_nopar`

```r
Omega <- Omega_est_nopar(data_list=summarydata,
                         ldscore.dir = "./eur_w_ld_chr")
```

The input from summary statistics:
  
- **summarydata**: a list of GWAS-summary-level data for exposure (*dat_x*) and two outcomes (*dat_y1*,*dat_y2*)
  
   *dat_x*: GWAS summary-level data for exposure, including

    1. rs number,
    2. effect allele,
    3. the other allele,
    4. sample size,
    5. a signed summary statistic (used to calculate z-score).

   For example, the *dat_x* with 3 SNPs can be represented as follows:
  
              SNP           b          se      frq_A1  A1 A2    P      N
       1: rs144155419 -0.0001011242 0.04969397  0.010  A  G 0.9983764 19563
       2:  rs58276399 -0.0040463990 0.01609338  0.111  C  T 0.8014823 19125
       3: rs141242758 -0.0062410460 0.01608838  0.111  C  T 0.6980775 19179

   *dat_yk*: GWAS summary-level data for the $k$-th outcome which is similar as *dat_x*.
 
- **ldscore.dir**: specify the path to the LD score files.
- **nCores**: The number of required cores or *NA*.
- **system_used**: The system used.

The argument **ldscore.dir** specifies the path to LD score files. Because the GWASs used for this example are based on European samples, we can use the LD score files from <https://github.com/yuanzhongshang/MAPLE/tree/main/example/eur_w_ld_chr>, which are provided by the ldsc software (<https://github.com/bulik/ldsc>). These LD Scores were computed using 1000 Genomes European data. Users can also calculate the LD scores by themselves.

Users can specify the rs number, effect allele, and the other allele using the arguments "snp_col," "A1_col," and "A2_col," respectively. Users may designate one or both of the following columns for calculating z-scores: "b_col" (effect size), "se_col" (standard error), "z_col" (z-score), and "p_col" (p-value). The sample size can be defined using the "n_col" argument. Alternatively, in the absence of a designated sample size column, users can utilize the "n" argument to indicate the total sample size for each SNP. Incorporating the minor allele frequency ("freq_col") column, if available, is advisable as it aids in filtering out low-quality SNPs.

The functions `Omega_est` and `Omega_est_nopar` will also conduct the following quality control procedures:
  
- extract SNPs in HapMap 3 list,
- remove SNPs with minor allele frequency $< 0.05$ (if freq_col column is available),
- remove SNPs with alleles not in (G, C, T, A),
- remove SNPs with ambiguous alleles (G/C or A/T) or other false alleles (A/A, T/T, G/G or C/C),
- exclude SNPs in the complex Major Histocompatibility Region (Chromosome 6, 26Mb-34Mb),
- remove SNPs with $\chi^2 > \chi^2_{max}$. The default value for $\chi^2_{max}$ is $max(N/1000, 80)$.

Now, we can check the estimates with the following commands:
  
```r
Omega
#$Omega
#            [,1]        [,2]        [,3]
#[1,]  1.06320820 -0.04611552 -0.01709125
#[2,] -0.04611552  1.02055474  0.49859715
#[3,] -0.01709125  0.49859715  1.03129944

#$Omega_se
#           [,1]       [,2]       [,3]
#[1,] 0.03158760 0.02169158 0.01934258
#[2,] 0.02169158 0.03412790 0.02355318
#[3,] 0.01934258 0.02355318 0.04228503
```

The output contains:
  
- **Omega**: the estimate of $\mathbf\Omega$, the off-diagonal elements of `Omega` are the intercept estimates of cross-trait LD score regression; the diagonal elements of `Omega` are the intercept estimates of single-trait LD score regressions.

- **Omega.se**: the estimated matrix consists of the standard errors of the intercept estimates obtained from LD score regression.

Users have the option to skip this step and set the estimate `Omega` of $\mathbf\Omega$ to the identity matrix if there is no confounding arising from sample structure.

### Step 2: Running METEOR

The `METEOR` function utilizes a scalable sampling-based algorithm to acquire calibrated $p$-values.

```r
#load the z-score
zscorex =fread(paste0("./zscorex.txt"),head=F)
zx<-as.matrix(zscorex,ncol=1)
zscorey_1 =fread(paste0("./zscorey1.txt"),head=F)
zy1<-as.vector(zscorey_1[[1]])
zscorey_2 =fread(paste0("./zscorey2.txt"),head=F)
zy2<-as.vector(zscorey_2[[1]])
Zscore <- cbind(zx,zy1,zy2)
Zscore <- t(Zscore)
#load the LD matrix
sigma<-fread(paste0("./Sigma.txt"),head=F)
Sigma <- as.matrix(sigma)
#load the sample size
N <- matrix(c(50000,50000,50000),ncol=1)
#load the correlation matrix and corresponding standard errors
load("./Omega.rda") #file name: Omega
Omega_est <- Omega$Omega
Omega_se <- Omega$Omega_se

result<-METEOR(Zscore,Sigma,N,Omega_est,Omega_se,Gibbsnumber=1000,burninproportion=0.2,
               pi_beta_shape=0.5,pi_beta_scale=4.5,pi_1_shape=0.5,pi_1_scale=1.5,
               pi_0_shape=0.05,pi_0_scale=9.95)
```

The input from summary statistics:
  
- **Zscore**: the Zscore with $K+1$ rows and p columns of the SNP effect size matrix for the exposure and $K$ outcomes.
- **Sigma**: the LD matrix for the SNPs selected from the exposure can be obtained by using the weighted average LD matrix from $K$ outcomes. If individual data is unavailable, the LD matrix can also be derived from a reference panel.
- **N**: the sample sizes of exposure and $K$ outcomes GWASs.
- **Omega_est**: the correlation matrix derived by LDSC.
- **Omega_se**: the standard errors of elements of $Omega$ derived by LDSC.
- **Gibbsnumber**: the number of Gibbs sampling iterations with the default to be 1000.
- **burninproportion**:  the proportion to burn in from Gibbs sampling iterations, with default to be 20%.
- **lambda**: the tuning parameter used to ensure that the correlation matrix is invertible.
- **pi_beta_shape**: the prior shape paramter for $\pi_\beta$ with the default to be 0.5.
- **pi_beta_scale**: the prior scale paramter for $\pi_\beta$ with the default to be 4.5.
- **pi_1_shape**: the prior shape paramter for $\pi_1$ with the default to be 0.5.
- **pi_1_scale**: the prior scale paramter for $\pi_1$ with the default to be 1.5.
- **pi_0_shape**: the prior shape paramter for $\pi_0$ with the default to be 0.05.
- **pi_0_scale**: the prior scale paramter for $\pi_0$ with the default to be 9.95.

Note that,  we use $p=5\times{10}^{-8}$ for METEOR to select candidate IVs without LD clumping. However, if the number of SNPs is too much, such as (greater than 10000), suggesting using LD clumping with $r^2=0.5$ to select candidate IVs. Additionally, users can employ the LD matrix derived from the weighted average LD matrix from k outcomes, or an LD reference panel as **Sigma**, provided that no additional LD matrices are available for the SNPs in the $k$ outcome data.

Now, we can check the estimates from METEOR:
  
```r
result$causal_effect
#          [,1]       [,2]
#[1,] 0.0785606 0.01103452

result$causal_pvalue_single
#             [,1]     [,2]
#[1,] 1.577082e-05 0.595196

result$causal_pvalue_overall
#             [,1]
#[1,] 6.608792e-05
```

The output from `METEOR` is a list containing:
  
- **causal_effect**: the estimate of $K$ causal effects.
- **causal_pvalue_single**: the $p$ values for the causal effects in single tests.
- **causal_pvalue_overall**: the $p$ value for the causal effect in overall tests.
- **cause_sd**: the standard deviation for the causal effects.
- **cause_cov**: the covariance for the causal effects.
- **sigmabeta**: the variance estimate for the SNP effect sizes on the exposure.
- **sigmaeta**: the variance estimates for the horizontal pleiotropy effects.
- **Omega**: the correlation matrix for one exposure and $K$ outcomes.
- **pi_beta**: the proportion of selected SNPs, which show non-zero effects on exposure.
- **pi_1**: the proportion of selected SNPs showing horizontal pleiotropy.
- **pi_0**: the proportion of non-selected SNPs showing horizontal pleiotropy.
