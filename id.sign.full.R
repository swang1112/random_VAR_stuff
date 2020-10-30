Rcpp::sourceCpp("sign_fast.cpp")
require(ggplot2)
require(dplyr)

#' Identification by Sign Restrictions
#'
#' @param Model A class of 'varest' or 'vec2var'
#' @param iter Integer, number of successful draws
#' @param Signmat Matrix of dimension (KxK), sign restrictions for impulse response functions:
#'                1 for positive; -1 for negative; 0 for unrestricted
#' @param r_ahead Integer, restricted horizons of impulse response functions: 
#'                0 if only impact effects are subject restrictions
#' @param n_ahead Integer, horizon of IRFs 
#' @param Ci Double, quantiles of accepted set value between [0,1] 
#' @param Epsname Char, shock names, by default Epsname = Yname
#' @param Plot Logica, Plot = TRUE for visualization of impulse response functions
#' @param Core Integer, number of cores available for parallel computing
#'
#' @return List
#' @export id.sign.full
#' @author ©Shu Wang 2020
#' 
#' @examples NULL
id.sign.full = function(Model, iter = 1000, Signmat, r_ahead, n_ahead = 15, 
                   Ci = 0.68, Epsname = NULL, Plot = TRUE, Core){
  
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
  SR_kernel(Accept_model, K = K, C = C, Signmat_r = Signmat, r_ahead = r_ahead, 
            A_hat = A_hat, iter = iter, n_ahead = n_ahead+1, Core = Core)
  
  Accept_model_flat = matrix(0, nrow = iter, ncol = K^2*(n_ahead+1)) 
  for (h in 1:(n_ahead+1)) {
    Accept_model_flat[,((K^2*(h-1)+1):(K^2*h))] = Accept_model %>% lapply("[[", h) %>% unlist %>% matrix(nrow = iter, ncol = K^2, byrow = T)
  }
  
  Medians = Accept_model_flat %>% apply(2, median)
  Sd      = Accept_model_flat %>% apply(2, sd)
  P_quant = Accept_model_flat %>% apply(2, function(x){quantile(x,probs = c((1-Ci)/2, (1+Ci)/2))})
  
  # standardzation for median target (MT) selection
  ## avoid dividing zero
  Sd[which( Sd < 1e-8 )] = 1
  Accept_model_stand = Accept_model_flat %>% apply(1, function(x){ (x - Medians)/Sd }) %>% t
  MT = which.min(Accept_model_stand[,1:(K^2*(r_ahead+1))] %>% apply(1, function(x){crossprod(x)}))
  
  erg = list()
  erg[["B"]] = Accept_model[[MT]][[1]]
  erg[["A_hat"]] = A_hat
  erg[["eps"]] = solve(Accept_model[[MT]][[1]]) %*% t(u)
  erg[["IRFs"]] = Accept_model[[MT]]
  
  # IRFs Plot ---------------------------------------------------------------
  if(Plot){
    
    # convert things to datafram for ggplot
    ready2plot = function(x, value_name = "value"){
      
      Out           = matrix(0, nrow = n_ahead+1, ncol = K^2 + 1)
      colnames(Out) = rep("h", K^2 + 1)
      Out[,1]       = 0:n_ahead
      c             = 1
      temp          = array(x, dim = c(K, K, n_ahead+1))
      for (i in 1:K) {
        for (j in 1:K) {
          c = c + 1
          Out[,c] = temp[i,j,]
          colnames(Out)[c] = paste("epsilon[", Epsname[j],"]", "%->%", Yname[i])
        }
      }
      Out %>% as.data.frame %>% reshape2::melt(id = "h", value.name = value_name)
    }
    
    # prepare to plot
    Response.MT = ready2plot(Accept_model_flat[MT,])
    Response.M  = ready2plot(Medians)
    Response.L  = ready2plot(P_quant[1,], "L")
    Response.U  = ready2plot(P_quant[2,], "U")
    
    Response.Ci   = Response.L %>% left_join(Response.U, by = c("variable", "h"))
    Response.plot = Response.M %>% left_join(Response.Ci, by = c("variable", "h")) %>% tibble::add_column(Label = "median")
    Response.plot = rbind(Response.plot, Response.MT %>% left_join(Response.Ci, by = c("variable", "h")) %>% tibble::add_column(Label = "median target"))
    
    IRF_plot = Response.plot %>% ggplot(aes(x = h, y = value, group = Label)) +
      geom_line(aes(linetype = Label), alpha = 1, size = 0.6) +
      geom_hline(yintercept = 0, color = 'red') +
      facet_wrap(~variable, scales = "free_y", labeller = label_parsed, ncol = K) +
      geom_ribbon(aes(ymin = L, ymax = U), alpha = 0.26) +
      xlab(" ") + ylab(" ") +
      theme_bw() + theme(legend.position = "none")
    
    erg[["Response"]] = Response.plot
    plot(IRF_plot)
  }
  
  return(erg)

}


