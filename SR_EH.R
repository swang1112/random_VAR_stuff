SR_EH = function(Model, iter = 1000, num_slow = 2, target = 3, 
                 Signmat_0, Signmat_r, r_start = 1, r_end = 3, n_ahead = 15, 
                 Ci = 0.68, Epsname = NULL){
  
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
  SR_EH(Accept_model, K = K, num_slow = num_slow, C = C, target = target, 
        Signmat_0 = Signmat_0, Signmat_r = Smat, r_start = r_start, r_end = r_end, 
        iter = iter, A_hat = A_hat, n_ahead = n_ahead)

# -------------------------------------------------------------------------

  
  n_ahead = 5
  iter = 10
SR_EH(Accept_model, K = K, num_slow = 1, C = C, target = 3, Signmat_0 = Signmat_0, Signmat_r = Smat, 
      r_start = 1, r_end = 1, iter = 10, A_hat = A_hat, n_ahead = 5)

# -------------------------------------------------------------------------


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
  
  IRF_combi = data.frame( "MT" = Accept_model_flat[MT,], "Median" = Medians,
                         "L" = P_quant[1,], "U" = P_quant[2,])
  Labs = rep("h", K^2*n_ahead)
  h = rep(0, K^2*n_ahead)
  c = 1
  for (m in 1:n_ahead) {
    for (i in 1:K) {
      for (j in 1:K) {
        Labs[c] = paste("epsilon[", Epsname[i],"]", "%->%", Yname[j])
        h[c]    = m
        c       = c + 1
      }
    }
  }
  
  IRF_combi = cbind(h, Labs, IRF_combi) %>% reshape2::melt(id = c("h", "Labs"))
  
  IRF_combi$Labs = factor(IRF_combi$Labs, levels = rev.default(IRF_combi$Labs %>% levels))
  
  ggplot(data = IRF_combi %>% filter(variable == "MT"), aes(x = h, y = value, group = Labs)) + theme_bw() +
    geom_line( alpha = 0.9) +
    geom_hline(yintercept = 0, color = 'red') +
    facet_wrap(~Labs, scales = "free_y", labeller = label_parsed) 
    geom_ribbon(aes(ymin = L, ymax = U), alpha = 0.2) +
    xlab("Horizon") + ylab("Response") +
    theme_bw()
  
  
  erg <- list()
  erg[["B"]] = Accept_model[[MT]][[1]]
  erg[["A_hat"]] = A_hat
  erg[["eps"]] = solve(Accept_model[[MT]][[1]]) %*% t(u)
  erg[["IRFs"]] = Accept_model[[MT]]
} 


