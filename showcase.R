rm(list = ls())
source("SR_EH.R")

library(svars)
data(USA)
var1 <- VAR(USA, p = 3, "const")

# restrictions on impact-multiplier
Signmat_0 = matrix(c(1, 0, 0, 
                     1, 1, -1, 
                     -1, 1, 1), nrow = 3, ncol = 3, byrow = T)

# restrictions on IRFs (only the third)
Smat = matrix(c(NaN, NaN, -1, 
                NaN, NaN, -1, 
                NaN, NaN, NaN), nrow = 3, ncol = 3, byrow = T) 

# see "SR_EH.R" files for details!
result = SR_EH(Model = var1, 
               iter = 10000, 
               num_slow = 1, 
               target = 3, 
               Signmat_0 = Signmat_0,
               Signmat_r = Smat, 
               r_start = 1, 
               r_end = 1, 
               n_ahead = 15, 
               Ci = 0.8, Plot = FALSE)

# (Median Target) MT impact multiplier
result$B

# MT IRFs
result$IRFs

# MT Monetary policy shock
result$eps

Rcpp::sourceCpp("sign_EH.cpp")
Smat = matrix(c(NaN, NaN, NaN, 
                NaN, 1, -1, 
                NaN, NaN, 1), nrow = 3, ncol = 3, byrow = T) 
Signmat_0 = matrix(c(1, 0, 0, 
                     1, -3, -1, 
                     -1, 1, 1), nrow = 3, ncol = 3, byrow = T)

agno_sign_check(Signmat_0, Smat)

