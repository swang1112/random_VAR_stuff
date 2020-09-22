#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins("cpp11")]]

using namespace Rcpp;

// choose 2 int from N (as index)
arma::imat choose2_fast(int N)
{
  arma::imat Out(N*(N-1)/2, 2);
  int c = 0;
  for(int i = 0; i < (N-1); i++)
  {
    for(int j = i+1; j < N; j++)
    { 
      Out(c,0) = i;
      Out(c,1) = j;
      c++;
    }
  }
  return Out;
}

// Givens rotation
// [[Rcpp::export]]
arma::mat givensQ_fast(arma::vec thetas, arma::mat B, int K)
{
  arma::mat Out = arma::eye(K, K);
  arma::imat Cmat = choose2_fast(K);
  for(int i = 0; i < K*(K-1)/2; i++)
  {

    arma::mat temp = arma::eye(K, K);
    temp(Cmat(i,0), Cmat(i,0)) = cos(thetas(i));
    temp(Cmat(i,1), Cmat(i,1)) = cos(thetas(i));
    temp(Cmat(i,0), Cmat(i,1)) = - sin(thetas(i));
    temp(Cmat(i,1), Cmat(i,0)) = sin(thetas(i));
    
    Out = Out * temp;
  }
  return B * Out.t();
}




