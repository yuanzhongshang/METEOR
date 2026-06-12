//library('Rcpp')
#include <RcppDist.h>
// [[Rcpp::depends(RcppArmadillo, RcppDist)]]


#include <RcppArmadillo.h>

#include <R.h>
#include <Rmath.h>
#include <cmath>
#include <stdio.h>
#include <stdlib.h> 
#include <cstring>
#include <ctime>
#include <Rcpp.h>

// Enable C++11 via this plugin (Rcpp 0.10.3 or later)
// [[Rcpp::plugins(cpp11)]]

using namespace Rcpp;
using namespace arma;
using namespace std;


//*******************************************************************//
//                MAIN FUNC                        //
//*******************************************************************//
//' METEOR
//' 
//' @export
// [[Rcpp::export]]

SEXP METEOR_CPP(SEXP Zscore_in,SEXP Sigma_in,SEXP N_in,SEXP p_in,SEXP k_in,SEXP Omega_in,SEXP Gibbsnumberin,
              SEXP burninproportion,SEXP pi_beta_shape_in,SEXP pi_beta_scale_in,
              SEXP pi_1_shape_in,SEXP pi_1_scale_in,SEXP pi_0_shape_in,SEXP pi_0_scale_in){// *
try{
	const int Gibbs_number = Rcpp::as<int>(Gibbsnumberin);
	const double burnin_p = Rcpp::as<double>(burninproportion);
  const int p = Rcpp::as<int>(p_in);
  const int k = Rcpp::as<int>(k_in);
	const double lamda_beta1 = Rcpp::as<double>(pi_beta_shape_in); 
	const double lamda_beta2 = Rcpp::as<double>(pi_beta_scale_in); 
	const double lamda_21 = Rcpp::as<double>(pi_1_shape_in); 
	const double lamda_22 = Rcpp::as<double>(pi_1_scale_in); 
	const double lamda_31 = Rcpp::as<double>(pi_0_shape_in); 
	const double lamda_32 = Rcpp::as<double>(pi_0_scale_in); 
	const arma::vec N = as<arma::vec>(N_in);
  const arma::mat Zscore = as<arma::mat>(Zscore_in);
  const arma::mat Sigma = as<arma::mat>(Sigma_in);
  const arma::mat Omega=as<arma::mat>(Omega_in);

  double sigma2beta_prior_shape=p/10+1;
  double sigma2beta_prior_scale=0.2;
  double sigma2beta = 2.0/p;
  double sigma2eta_prior_shape=p/5.0+1;
  double sigma2eta_prior_scale=0.2;
  double pi_beta = 0.1;
  vec xi = sqrt(N-1);

  vec onek = ones<vec>(k);
  vec onep = ones<vec>(p);
  vec onek1 = ones<vec>(k+1);
  vec sigma2eta = 1.0/p*onek;
  vec pi_1 = 0.25*onek;
  vec pi_0 = 0.005*onek;
  vec alpha = zeros<vec>(k+1);
  alpha(0) = 1;
  mat inv_Omega = inv(Omega,inv_opts::allow_approx);
  
  double U_beta,U_eta;
  vec latent_gamma = zeros<vec>(p);
  vec latent_beta = zeros<vec>(p);  
  mat latent_tau = zeros<mat>(k,p);
  mat latent_eta = zeros<mat>(k+1,p);
  
  double part_1,part_2,Probability_gamma,randu_number;
  double part_1_pleiotropy,part_2_pleiotropy,Probability_gamma_pleiotropy,randu_number_pleiotropy;
  int burnin = ceil(burnin_p*Gibbs_number);
  vec sample_pi_beta(Gibbs_number);
  mat sample_pi_1(k,Gibbs_number);
  mat sample_pi_0(k,Gibbs_number);
  mat sample_alpha(k,Gibbs_number);
  vec sample_sigma2beta(Gibbs_number);
  mat sample_sigma2eta(k,Gibbs_number);
  double post_sigma2beta_scale;
  double post_sigma2beta_shape;
  double post_sigma2eta_scale;
  double post_sigma2eta_shape;
  mat nn = kron(onep,xi.t());
  
  for(int m=0; m<Gibbs_number; ++m){
    
    for (int i=0; i<k; ++i){
      
      double sample_var_eta = 1.0/as_scalar((N(i+1)-1)*inv_Omega(i+1,i+1) + 1/sigma2eta(i));
     
      for(int j=0; j<p; ++j){
        
        double sample_var_beta =1.0/(1.0/sigma2beta+as_scalar((xi%alpha).t()*inv_Omega*(xi%alpha)));
        U_beta = as_scalar(((Zscore.col(j)).t()-Sigma.col(j).t()*latent_beta*((xi%alpha).t())-Sigma.col(j).t()*(latent_eta.t()%nn)
                              +latent_beta(j)*((xi%alpha).t()))*inv_Omega*(xi%alpha)*sample_var_beta);
        part_1 = exp(0.5*U_beta*U_beta/sample_var_beta+0.5*log(sample_var_beta)-0.5*log(sigma2beta)+log(pi_beta)+ 
          sum(latent_tau.col(j)%log(pi_1))+sum((1-latent_tau.col(j))%log(1-pi_1)));
       
        part_2 = exp(log(1-pi_beta)+ sum(latent_tau.col(j)%log(pi_0))+sum((1-latent_tau.col(j))%log(1-pi_0)));
        Probability_gamma = part_1/(part_1+part_2);
        randu_number = as_scalar(randu(1));
        
        if(randu_number <= Probability_gamma){
          latent_gamma(j) = 1;
          latent_beta(j) = as_scalar(randn(1)*sqrt(sample_var_beta)+U_beta);
        } else {
          latent_gamma(j) = 0;
          latent_beta(j) = 0;
        }
        
        U_eta = as_scalar(xi(i+1)*(Zscore.col(j).t()-latent_beta.t()*Sigma.col(j)*((xi%alpha).t())-
          Sigma.col(j).t()*(latent_eta.t()%nn))*inv_Omega.col(i+1)+
          (N(i+1)-1)*latent_eta(i+1,j)*inv_Omega(i+1,i+1))*sample_var_eta;
        part_1_pleiotropy = exp(0.5*U_eta*U_eta/sample_var_eta + 0.5*log(sample_var_eta)-0.5*log(sigma2eta(i))+
          latent_gamma(j)*log(pi_1(i))+ (1-latent_gamma(j))*log(pi_0(i)));
        part_2_pleiotropy = exp(latent_gamma(j)*log(1-pi_1(i))+ (1-latent_gamma(j))*log(1-pi_0(i)));
        Probability_gamma_pleiotropy = part_1_pleiotropy/(part_1_pleiotropy+part_2_pleiotropy);
        randu_number_pleiotropy = as_scalar(randu(1));
        if(randu_number_pleiotropy<= Probability_gamma_pleiotropy){
          latent_tau(i,j) = 1;
          latent_eta(i+1,j) = as_scalar(randn(1)*sqrt(sample_var_eta)+ U_eta);
        } else {
          latent_tau(i,j) = 0;
          latent_eta(i+1,j) = 0;
        }		
  
      }
      
      vec T2 = Sigma*latent_beta;
      double M1 = as_scalar(T2.t()*latent_beta);
      double indicator_alpha=as_scalar((N(i+1)-1)*inv_Omega(i+1,i+1))*M1;
      
      if(indicator_alpha > 0){
        double sample_alpha_var = 1.0/indicator_alpha;
        double sample_alpha_mean = as_scalar(xi(i+1)*(latent_beta.t()*Zscore.t()-T2.t()*(latent_eta.t()%nn)
                                              -M1*(xi%alpha).t())*inv_Omega.col(i+1)*sample_alpha_var+alpha(i+1));
        sample_alpha(i,m) = as_scalar(randn(1)*sqrt(sample_alpha_var)+ sample_alpha_mean);
        alpha(i+1)=sample_alpha(i,m);
      } else {
        sample_alpha(i,m) = 0;
        alpha(i+1)=sample_alpha(i,m);
      }
      
      double shape1_pi_1 = sum(latent_gamma % (latent_tau.row(i).t()))+lamda_21;
      double shape2_pi_1 = sum(latent_gamma % (1.0-latent_tau.row(i).t()))+lamda_22;
      sample_pi_1(i,m)= r_4beta(shape1_pi_1,shape2_pi_1,0,1);
      pi_1(i) = sample_pi_1(i,m);
      
      double shape1_pi_0 = sum((1-latent_gamma) % (latent_tau.row(i).t()))+lamda_31;
      double shape2_pi_0 = sum((1-latent_gamma) % (1.0-(latent_tau.row(i).t())))+lamda_32;
      sample_pi_0(i,m)= r_4beta(shape1_pi_0,shape2_pi_0,0,1);
      pi_0(i) = sample_pi_0(i,m);
      
      
      post_sigma2eta_scale = 2/(as_scalar(sum((latent_tau.row(i).t())%(latent_eta.row(i+1).t())%(latent_eta.row(i+1).t()))+2*sigma2eta_prior_scale));
      post_sigma2eta_shape = 0.5*sum(latent_tau.row(i).t())+ sigma2eta_prior_shape;
      post_sigma2eta_scale = abs(post_sigma2eta_scale);
      sample_sigma2eta(i,m) = 1.0/as_scalar(randg(1, distr_param(post_sigma2eta_shape,post_sigma2eta_scale)));
      sigma2eta(i) = sample_sigma2eta(i,m);
    }
    double shape1 = sum(latent_gamma)+lamda_beta1;
    double shape2 = sum(1.0-latent_gamma)+lamda_beta2;
    sample_pi_beta(m)= r_4beta(shape1,shape2,0,1);
    pi_beta = sample_pi_beta(m);
  
    post_sigma2beta_scale = 2/(as_scalar(sum(latent_gamma % latent_beta % latent_beta)+2*sigma2beta_prior_scale));
    post_sigma2beta_shape = 0.5*sum(latent_gamma)+ sigma2beta_prior_shape;
    post_sigma2beta_scale = abs(post_sigma2beta_scale);
    sample_sigma2beta(m) = 1.0/as_scalar(randg(1, distr_param(post_sigma2beta_shape,post_sigma2beta_scale)));
    sigma2beta = sample_sigma2beta(m);
  }
  
  double sigma2beta_estimate = mean(sample_sigma2beta.subvec(burnin,(Gibbs_number-1)));
  double pi_beta_estimate = mean(sample_pi_beta.subvec(burnin,(Gibbs_number-1)));
  rowvec sigma2eta_estimate = mean(sample_sigma2eta.submat(0,burnin,(k-1),(Gibbs_number-1)).t());
  rowvec pi_1_estimate = mean(sample_pi_1.submat(0,burnin,(k-1),(Gibbs_number-1)).t());
  rowvec pi_0_estimate = mean(sample_pi_0.submat(0,burnin,(k-1),(Gibbs_number-1)).t());
  rowvec alpha_estimate = mean(sample_alpha.submat(0,burnin,(k-1),(Gibbs_number-1)).t());
  rowvec alpha_sd_estimate = stddev(sample_alpha.submat(0,burnin,(k-1),(Gibbs_number-1)).t());
  mat alpha_cov_estimate = cov(sample_alpha.submat(0,burnin,(k-1),(Gibbs_number-1)).t());

  return List::create(Rcpp::Named("alpha") = alpha_estimate,
                      Rcpp::Named("alpha_sd") = alpha_sd_estimate,
                      Rcpp::Named("alpha_cov") = alpha_cov_estimate,
                      Rcpp::Named("sigma2beta") = sigma2beta_estimate,
                      Rcpp::Named("sigma2eta") = sigma2eta_estimate,
                      Rcpp::Named("pi_beta") = pi_beta_estimate,
                      Rcpp::Named("pi_1") = pi_1_estimate,
                      Rcpp::Named("pi_0") = pi_0_estimate);
                             
	} catch( std::exception &ex ) {
		forward_exception_to_r( ex );
	} catch(...) {
		::Rf_error( "C++ exception (unknown reason)..." );
	}
	return R_NilValue;
}// end func
