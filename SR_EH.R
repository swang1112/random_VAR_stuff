require(dplyr)
require(ggplot2)
Rcpp::sourceCpp("sign_EH.cpp")

#' Eickmeier und Hofmann (2013) 's Scheme of Cholesky and Sign Restrictions
#'
#' @param Model A class of 'varest' or 'vec2var'
#' @param iter Integer, number of successful draws
#' @param num_slow Integer, number of slow-moving variables, which MP shock does not instantaneously affects
#' @param target Integer, (column) position of MP shock 
#' @param Signmat_0 Matrix of dimension (KxK), sign restrictions for impact multiplier. Achtung! This needs to be a lower-triangular matrix!!
#' @param Signmat_r Matrix of dimension (KxK), sign restrictions for IRFs
#' @param r_start Integer, starting point of restricted IRFs horizon, by default r_start = 1
#' @param r_end Integer, ending point of restricted IRFs horizon
#' @param n_ahead Integer, horizon of IRFs 
#' @param Ci Double, confidence interval value between [0,1] (bootstrap under development)
#' @param Epsname Char, shock names, by default Epsname = Yname
#' @param Plot 
#' @param Core Integer, Number of cores available for parallel computing
#'
#' @return List
#' @export SR_EH
#'
#' @examples NA
SR_EH = function(Model, iter = 1000, num_slow = 2, target = 3, 
                 Signmat_0, Signmat_r, r_start = 1, r_end = 3, n_ahead = 15, 
                 Ci = 0.68, Epsname = NULL, Plot = TRUE, Core = 8){
  
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
  SR_EH_kernel(Accept_model, K = K, num_slow = num_slow, C = C, target = target, 
        Signmat_0 = Signmat_0, Signmat_r = Smat, r_start = r_start, r_end = r_end, 
        iter = iter, A_hat = A_hat, n_ahead = n_ahead, core = Core)
  
  Accept_model_flat = matrix(0, nrow = iter, ncol = K^2*n_ahead) 
  for (h in 1:n_ahead) {
    Accept_model_flat[,((K^2*(h-1)+1):(K^2*h))] = Accept_model %>% lapply("[[", h) %>% unlist %>% matrix(nrow = iter, ncol = K^2, byrow = T)
  }
  
  Medians = Accept_model_flat %>% apply(2, median)
  Sd      = Accept_model_flat %>% apply(2, sd)
  P_quant = Accept_model_flat %>% apply(2, function(x){quantile(x,probs = c((1-Ci)/2, (1+Ci)/2))})
  
  # standardzation for median target (MT) selection
  ## avoid dividing zero
  Sd[which( Sd < 1e-8 )] = 1
  Accept_model_stand = Accept_model_flat %>% apply(1, function(x){ (x - Medians)/Sd }) %>% t
  MT = which.min(Accept_model_stand[,(1:(K^2*r_end))] %>% apply(1, function(x){crossprod(x)}))

  erg = list()
  erg[["B"]] = Accept_model[[MT]][[1]]
  erg[["A_hat"]] = A_hat
  erg[["eps"]] = solve(Accept_model[[MT]][[1]]) %*% t(u)
  erg[["eps"]] = erg[["eps"]][target,]
  erg[["IRFs"]] = Accept_model[[MT]]

  
  # IRFs Plot ---------------------------------------------------------------
  if(Plot){

    # convert things to datafram for ggplot
    ready2plot = function(x, value_name = "value"){
      
      Out           = matrix(0, nrow = n_ahead, ncol = K^2 + 1)
      colnames(Out) = rep("h", K^2 + 1)
      Out[,1]       = 1:n_ahead
      c             = 1
      temp          = array(x, dim = c(K, K, n_ahead))
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
    
    # Partial identification
    Response.partial <- data.frame()
    for (k in 1:(2*K)) {
      Response.temp <- Response.plot[((n_ahead)*K*(k-1) + 1): ((n_ahead)*K*k),]
      Response.partial <- bind_rows(Response.partial,Response.temp[((n_ahead)*(target-1)+1):((n_ahead)*target),])
    }
    
    IRF_plot = Response.partial %>% ggplot(aes(x = h, y = value, group = Label)) +
      geom_line(aes(linetype = Label, color = Label), alpha = 0.9) +
      geom_hline(yintercept = 0, color = 'red') +
      facet_wrap(~variable, scales = "free_y", labeller = label_parsed, ncol = 1) +
      geom_ribbon(aes(ymin = L, ymax = U), alpha = 0.2) +
      xlab("Horizon") + ylab("Response") +
      theme_bw() 
    
    erg[["Response"]] = Response.plot
    plot(IRF_plot)
  }
  
  return(erg)
} 


