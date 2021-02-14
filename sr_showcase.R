rm(list = ls())
source('id.sr_stand.R')
library(svars)
data(USA)

USA %>% head
varm = VAR(USA, p = 2, type = "const")
Signmat = matrix(c(1,-1, NA,
                   1, 1, -1,
                   1, 1,  1), nrow = 3, ncol = 3, byrow = T)

svar_sr = id.sr(varm, accept = 1000, Signmat = Signmat, n_ahead = 20, Core = parallel::detectCores()-1)
