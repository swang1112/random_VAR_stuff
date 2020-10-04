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
Smat = matrix(c(NA, NA, -1, 
                NA, NA, -1, 
                NA, NA, 1), nrow = 3, ncol = 3, byrow = T) 

# see "SR_EH.R" files for details!
result = SR_EH(Model = var1, 
               iter = 1000, 
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

