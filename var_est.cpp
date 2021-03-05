#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]

/*
// [[Rcpp::export]]

arma::mat var_est(arma::mat y, arma::mat Z, int& K){
  
  arma::mat ZZ = arma::inv(Z*arma::trans(Z)) * Z; 
  arma::mat I  = arma::eye(K, K);
  arma::mat Out = arma::kron(ZZ, I) * y;
  return Out;
  
}

// [[Rcpp::export]]

arma::mat var_resid(arma::mat y, arma::mat Z, arma::mat beta, int& K){
  
  arma::mat I  = arma::eye(K, K);
  arma::mat Out = y - arma::kron(arma::trans(Z), I) * beta;
  return Out;
  
}
*/
 
// [[Rcpp::export]]

arma::mat var_est1(arma::mat y, arma::mat Z, int& K){
  
  arma::mat ZZ = arma::inv(Z*arma::trans(Z)) * Z; 
  arma::mat I  = arma::eye(K, K);
  arma::mat Out = arma::kron(I, ZZ) * y;
  return Out;
  
}

// [[Rcpp::export]]

arma::mat var_resid1(arma::mat y, arma::mat Z, arma::mat beta, int& K){
  
  arma::mat I  = arma::eye(K, K);
  arma::mat Out = y - arma::kron(I, arma::trans(Z)) * beta;
  return Out;
  
}


// [[Rcpp::export]]

arma::mat var_hess(arma::mat Z, arma::mat Cov_u){
  
  arma::mat Out = arma::kron(Cov_u, arma::inv(Z*arma::trans(Z)));
  return Out;
  
}

/*
// [[Rcpp::export]]
arma::mat var_sur(arma::mat y, arma::mat Z, int K){
  
  int T = y.n_rows / K;
  arma::mat ZZ = arma::inv(Z*arma::trans(Z)) * Z; 
  arma::mat Out = arma::zeros(Z.n_rows*K, 1);
  for (int i = 0; i < K; i ++)
  {
    Out.submat(Z.n_rows*i, 0, Z.n_rows*(i+1)-1, 0) =  arma::inv(Z*arma::trans(Z)) * Z * y.submat(T*i, 0, T*(i+1)-1, 0);
  }
  
  return Out;
  
}
*/