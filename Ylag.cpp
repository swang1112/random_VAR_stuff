#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]


// [[Rcpp::export]]

arma::mat Ylag(arma::mat y, int p, int c ){
  
  arma::mat YLag = arma::zeros(y.n_rows, y.n_cols * p);
  
  for(int i = 0; i < p; ++i){
    YLag.submat(i, i * y.n_cols, y.n_rows-1, i * y.n_cols + (y.n_cols-1)) = y.submat(0, 0, y.n_rows - i - 1, y.n_cols-1 );
  }
  
  arma::mat YLagOut = YLag.submat(p-1,0, YLag.n_rows -2, YLag.n_cols - 1);
  
  if (c == 0)
  {
    return YLagOut;
  } else if (c == 1)
  {
    arma::vec c1(YLagOut.n_rows, arma::fill::ones);
    arma::mat YLagOut1 = arma::join_horiz(c1, YLagOut);
    return YLagOut1;
  } else if (c == 2)
  {
    arma::vec c1(YLagOut.n_rows, arma::fill::ones);
    arma::vec c2 = arma::linspace(1, YLagOut.n_rows, YLagOut.n_rows);
    arma::mat YLagOut2 = arma::join_horiz(c1, c2, YLagOut);
    return YLagOut2;
  } else
  {
    arma::mat place_holder = arma::zeros<arma::mat>(1,1);
    return place_holder;
  }
  
}
