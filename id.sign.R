Rcpp::sourceCpp("sign_fast.cpp")

id.sign = function(Model, iter = 1000, Partial = NULL, Signmat_r, r_ahead, n_ahead = 15, 
                   Ci = 0.68, Epsname = NULL, Plot = TRUE){
  p     = Model$p 
  u     = Model %>% resid # residuals
  Tob   = u %>% nrow # observations
  K     = u %>% ncol # system dimension
  Covmt = crossprod(u) / (Tob - K*p - 1)
  C     = Covmt %>% chol %>% t # low-trig factor
  A_hat = vars::Bcoef(Model)[,1:(K*p)] # AR coeff
  Yname = names(Model$varresult) # Variable names
  if(is.null(Epsname)){Epsname = Yname}
  
  Accept_model = vector("list", length = iter)
  SR_kernel(Accept_model, K = K, num_slow = num_slow, C = C, target = target, 
               Signmat_0 = Signmat_0, Signmat_r = Smat, r_start = r_start, r_end = r_end, 
               iter = iter, A_hat = A_hat, n_ahead = n_ahead)
  
  
  
  
}