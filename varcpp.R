Rcpp::sourceCpp('funs/Ylag.cpp')
Rcpp::sourceCpp('funs/var_est.cpp')

#' LS estimation of VAR model
#'
#' @param Y matrix of dimension (TxK) with T observations and K endogeneous variables
#' @param p integer for the lag order
#' @param c integer for the deterministic term: c = 0 for no deterministic term, c = 1 for constant intercept, c = 2 for deterministic trend
#'
#' @return list
#' @export varcpp
#' @author Shu Wang \email{shu.wang@uni-goettingen.de}
#' 
#' @examples NULL
varcpp = function(Y, p, c = 1){
  K     = ncol(Y)
  totob = nrow(Y) 
  Z     = t(Ylag(Y, p, c))
  y     = matrix(Y[-c(1:p),], ncol = 1)
  b     = var_est1(y, Z, K)
  u     = var_resid1(y, Z, b, K)
  umat  = matrix(u, nrow = totob-p, ncol = K)
  Cov_u = crossprod(umat)/(totob-p-K*p-1)
  b_hess = var_hess(Z, Cov_u)
  #b_se  = sqrt(diag(b_hess))
  erg   = list('dat' = Y,
               'Z' = Z,
               'b' = b,
               'hess' = b_hess,
               'u' = umat,
               'Cov_u' = Cov_u,
               'A_hat' = matrix(b, nrow = K, byrow = T))
}
