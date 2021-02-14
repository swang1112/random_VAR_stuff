Rcpp::sourceCpp("IRF.cpp")
require(magrittr)
require(ggplot2)
#' fast standard sign restriction 
#' ©Shu Wang 2020
#'
#' @param Model Object of class 'varest' or 'vec2var'
#' @param accept Int, number of candiates
#' @param Signmat Mat, sign restrictions (1 for positive, -1 for negative, NA for agnostic/unrestricted)
#' @param n_ahead Int, horizon of impulse responses to be calculated
#' @param Ci Num, point-wise quantiles of accepted draws
#' @param Core Int, number of CPU kernels for parrallel computation
#' @param Plot Bool, plot IRFs if Plot = TRUE
#'
#' @return list
#' @export id.sr
#'
#' @author Shu Wang \email{shu.wang@uni-goettingen.de}
#'
#' @examples NULL
id.sr = function(Model, accept = 1000, Signmat, n_ahead = 30, Ci = 0.9, Core = 1, Plot = T){

  p     = Model$p 
  u     = Model %>% resid # residuals
  Tob   = u %>% nrow # observations
  K     = u %>% ncol # system dimension
  Covmt = crossprod(u) / (Tob - K*p - 1)
  C     = Covmt %>% chol %>% t # low-trig factor
  A_hat = vars::Bcoef(Model)[,1:(K*p)] # AR coeff
  Yname = names(Model$varresult) # Variable names
  
  Accept_model = vector("list", length = accept)
  Accept_model = pbapply::pblapply(Accept_model, .sr_stand, C = C, sr = Signmat, K = K, cl = Core)
  
  foo =  parallel::mclapply(Accept_model, .transformer, Partial = Partial, A_hat = A_hat, n_ahead = n_ahead+1, mc.cores = Core)
  new_dim = dim(foo[[1]])
  Accept_model_flat = matrix(0, nrow = accept, ncol = new_dim[1]*new_dim[2]) 
  for (i in 1:accept) {
    Accept_model_flat[i,] = foo[[i]] %>% t %>% matrix(nrow = 1)
  }
  
  Medians = Accept_model_flat %>% apply(2, median)
  Sd      = Accept_model_flat %>% apply(2, sd)
  P_quant = Accept_model_flat %>% apply(2, function(x){quantile(x,probs = c((1-Ci)/2, (1+Ci)/2))})
  
  Sd[which( Sd < 1e-8 )] = 1
  Accept_model_stand = Accept_model_flat %>% apply(1, function(x){ (x - Medians)/Sd }) %>% t

  MT = which.min(Accept_model_stand[,1:new_dim[2]] %>% apply(1, function(x){crossprod(x)}))
  
  erg = list()
  erg[['Accept_model']] = Accept_model
  erg[["B"]] = Accept_model[[MT]]
  erg[["Bcoef"]] = A_hat
  erg[["IRFs"]] = foo[[MT]]
  erg[["CB"]] = P_quant
  erg[["Sd"]] = Sd
  erg[["Median"]] = Medians
  erg[["Signmat"]] = Signmat
  
  if (Plot) {
    point_est = data.frame(foo[[MT]])
    c = 1
    for (j in 1:K) {
      for (i in 1:K) {
        colnames(point_est)[c] = paste("epsilon[", Yname[j],"]", "%->%", Yname[i])
        c = c + 1
      }
    }
    L = P_quant[1,] %>% matrix(ncol = ncol(point_est), byrow = T)
    U = P_quant[2,] %>% matrix(ncol = ncol(point_est), byrow = T)
    colnames(L) = colnames(U) = colnames(point_est)
    h = 0:(nrow(point_est)-1)
    L = cbind(h, L) %>% as.data.frame %>% reshape2::melt(id.vars = 'h', value.name = "L")
    U = cbind(h, U) %>% as.data.frame %>% reshape2::melt(id.vars = 'h', value.name = "U")
    
    IRF_plot = cbind(h, point_est) %>% as.data.frame %>% reshape2::melt(id.vars = 'h') %>% 
      dplyr::left_join(L, by = c('h', 'variable')) %>% 
      dplyr::left_join(U, by = c('h', 'variable')) %>%
      ggplot(aes(x = h, y = value)) + geom_line() + geom_hline(yintercept = 0, color = 'red') +
      facet_wrap(~variable, scales = "free_y", labeller =label_parsed, nrow = K, dir = "v") +
      geom_ribbon(aes(ymin = L, ymax = U), alpha=0.3) + xlab(" ") + ylab(" ") +
      theme_bw()
    
    erg[["IRF_plot"]] = IRF_plot
    show(IRF_plot)
  }
  
  return(erg)
}


# internal functions ------------------------------------------------------
.sr_stand = function(x, C, sr, K){
  pass = 0
  while (pass == 0) {
    G = matrix(rnorm(K*K), nrow = K, ncol = K)
    Q = qr.Q(qr(G))
    B = C %*% Q
    B = .positiv_diag(B)
    pass = .sign_check(B, sr)
  }
  B
}

.sign_check = function(x, sr){
  sign_test  <- sign(x) == sr
  all(sign_test, na.rm = TRUE)
}

.positiv_diag = function(x){
  for (kk in 1:ncol(x)) {
    if (x[kk,kk]<0) {
      x[,kk] = x[,kk] * -1
    }
  }
  x
}

.transformer = function(B, Partial, A_hat, n_ahead){
  IRF_fast(A_hat = A_hat, B_hat = B, horizon = n_ahead) %>% sapply(function(x) matrix(x, nrow = nrow(B), ncol = ncol(B))) %>% t() 
}

