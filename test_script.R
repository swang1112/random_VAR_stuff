
load(file = "data_K6.RData")

var_est1      <- vars::VAR(ds_lkp_ts[,-1], type = "both", p = 3)  #c(2,3,4,7,5,6)
colnames(ds_lkp_ts)
###______ Sign restrictions with zero restrictions:____________

SR_mat_list <- list()

### without money 
Signmat_0 <- matrix(NaN, nrow = 6, ncol = 6)
Signmat_0[c(1,2),c(3:6)] <- 0
Signmat_0[c(6),c(6)] <- 1
Signmat_0
SR_mat_list$K6_baseline1$Signmat_0 <- Signmat_0


#########################################################

### K6 without money 

Smat <- matrix(NaN, nrow = 6, ncol = 6)
Smat[c(6),6]     <-  1
Smat[c(1,2),6]   <- -1
Smat
SR_mat_list$K6_baseline1$Smat <- Smat


Smat <- matrix(NaN, nrow = 6, ncol = 6)
Smat[c(6),6]     <-  1
Smat[c(1,2),6]   <- -1
Smat
SR_mat_list$K6_baseline1$Smat <- Smat
Smat <- matrix(NaN, nrow = 6, ncol = 6)
Smat[c(6),6] <- 1
Smat[c(1,2,4,5),6] <- -1
Smat
SR_mat_list$K6_baseline_full$Smat <- Smat
Smat <- matrix(NaN, nrow = 6, ncol = 6)
Smat[c(6),6] <- 1
Smat[c(2),6] <- -1
Smat
SR_mat_list$K6_baseline_agnostic$Smat <- Smat


source("SR_EH.R")

Smat
Signmat_0


X = matrix(c(1,2,3,4,5,1, 
             1,2,3,4,-5,-4,
             rnorm(23),1), nrow = 6, ncol = 6, byrow = T)

X
agno_sign_check(X, Smat)


sr_results <- list()

r_start = 1
r_end = 5
Ci = 1

sr_results$agnostic_h1 <- SR_EH(Model = var_est1, 
                                
                                iter = 50, 
                                
                                num_slow = 2, 
                                
                                target = 6, 
                                
                                Signmat_0 = SR_mat_list$K6_baseline1$Signmat_0,
                                
                                Signmat_r = SR_mat_list$K6_baseline1$Smat, 
                                
                                r_start = r_start, 
                                
                                r_end = r_end, 
                                
                                n_ahead = 50, 
                                
                                Ci = Ci, Plot = TRUE, Core = 16)

for (h in (r_start+1):(r_end+1)) {
  print(sr_results$agnostic_h1$IRFs[[h]][1:2,6]<0)
}

