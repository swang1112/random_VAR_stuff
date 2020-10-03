Rcpp::sourceCpp("CvM_fast.cpp")



K = 3
B = matrix(rnorm(K^2), nrow = K, ncol = K)
theta = runif(choose(K,2), 0, pi)
all.equal(givens.Q.svars(theta, B), givensQ_fast(theta, B, K))




library(svars)
data(USA)
USA
p = 3
var1 <- VAR(USA, p = p, "const")

A = Bcoef(var1)[,1:(p*3)]
B = id.dc(var1)$B

B = diag(3)

B = chol(crossprod(resid(var1))/(174))

var1 %>% irf(n.ahead = 2, boot = F)
IrF(A, B, 2)

irf(id.dc(var1) , n.ahead = 2)$irf
array(dim = c(3,3,2))


"%^%" <- function(A, n){
  if(n == 0){
    diag(ncol(A))
  }else if(n == 1){
    A
  }else{
    A %*% (A %^% (n-1))
  }
}
# function to calculate impulse response
IrF <- function(A_hat, B_hat, horizon){
  k <- nrow(A_hat)
  p <- ncol(A_hat)/k
  if(p == 1){
    irfa <- array(0, c(k, k, horizon))
    for(i in 1:horizon){
      irfa[,,i] <- (A_hat%^%(i-1))%*%B_hat
    }
    return(irfa)
  }else{
    irfa <- array(0, c(k, k, horizon))
    irfa[,,1] <- B_hat
    Mm <- matrix(0, nrow = k*p, ncol = k*p)
    Mm[1:k, 1:(k*p)] <- A_hat
    Mm[(k+1):(k*p), 1 : ((p-1)*k)] <- diag(k*(p-1))
    Mm1 <- diag(k*p)
    for(i in 1:(horizon-1)){
      Mm1 <- Mm1%*%Mm
      irfa[,,(i+1)] <- Mm1[1:k, 1:k]%*%B_hat
    }
    return(irfa)
  }
}



Rcpp::sourceCpp("IRF.cpp")

n = 10
K = 5
B = matrix(rnorm(K^2), nrow = K, ncol = K)
all.equal(B%^%n, matexp(B, n))

eva = microbenchmark::microbenchmark("R" = B%^%n,
                                     "C++" = matexp(B, n),
                                     times = 1000)


eva
eva %>% autoplot()


library(svars)
library(ggplot2)
data(USA)
USA %>% head
p = 3
var1 <- VAR(USA, p = p, "const")
A_hat = Bcoef(var1)[,1:(3*p)]
B_hat = id.dc(var1)$B
IrF(A_hat, B_hat, 10)
irf(id.dc(var1), n.ahead = 3)
IRF_fast(A_hat, B_hat, 10)

eva = microbenchmark::microbenchmark("R" = IrF(A_hat, B_hat, 20),
                                     "Cpp" = IRF_fast(A_hat, B_hat, 20),
                                     times = 10000)
eva
autoplot(eva)

Rcpp::sourceCpp("sign_fast.cpp")

set.seed(1234)
A = matrix(c(0,rnorm(8)), nrow = 3, ncol = 3)
Smat  = matrix(c(NA, -1, 1, 1, -1, 1, 1, -1, -1), nrow = 3, ncol = 3)
sign_check(A, Smat) 
all(sign(A) == Smat)
eva = microbenchmark::microbenchmark("R" = all(sign(A) == Smat),
                                     "cpp" = sign_check(A, Smat), 
                                     times = 10000)
eva %>% autoplot()

typeof(out)
class(out)
out = matrix(c(1,1,0, -1, -1, 1, 0, 1, -1), nrow = 3, ncol = 3)

B = chol(A %*% t(A)) %>% t
B
solve(B)


Rcpp::sourceCpp("sign_fast.cpp")
library(svars)
library(dplyr)
data(USA)
USA %>% head
p = 3
var1 <- VAR(USA, p = p, "const")
A_hat = Bcoef(var1)[,1:(3*p)]
B_hat = id.dc(var1)$B
Smat = matrix(c(rep(NA, 6), -1, 1, 1), 3,3)
IRF_temp = IRF_fast(A_hat, B_hat, 10)
sign_check(IRF = IRF_temp, r_start = 2, r_end = 6, num_slow = 2, K = 3, Smat = Smat)
test_fast(3)




Rcpp::sourceCpp("sign_fast.cpp")
library(svars)
library(dplyr)
data(USA)
USA %>% head
p = 3
var1 <- VAR(USA, p = p, "const")
A_hat = Bcoef(var1)[,1:(3*p)]
B_hat = id.dc(var1)$B



Smat = matrix(c(rep(NA, 6), -1, 1, -1), 3,3)

Signmat_0 = matrix(c(-1, 1, 1, 1, -1, -1, 1, 1, -1), nrow = 3, ncol = 3)

var1 = VAR(USA, p = p, "const")
Model = var1
p     = Model$p
u     = Model %>% resid # residuals
Tob   = u %>% nrow # observations
K     = u %>% ncol # system dimension
Covmt = crossprod(u) / (Tob - K*p - 1)
C     = Covmt %>% chol %>% t # chol factor
A_hat = vars::Bcoef(Model)[,1:(K*p)] # reduced form parameters



num_slow = 1
thetas = runif((K-num_slow)*(K-num_slow-1)/2)
thetas = thetas * pi
R = givensQ_fast(thetas, K-num_slow)
Rotmat = diag(3)
Rotmat[((num_slow+1) :  K), ((num_slow+1) :  K)] = R
Rotmat
Bmat = solve(Rotmat %*% C)

for (j in 1:K) {
  if (Bmat[j,j] < 0) {
    Bmat[,j] <- Bmat[,j] * -1
  }
}

irf_temp = IRF_fast(A_hat, Bmat, 5)



Smat = matrix(c(-1, 1, 1, 
                1, 1, 1, 
                1, 1, 1), nrow = 3, ncol = 3, byrow = T)

Signmat_0 = matrix(c(1, 0, 0, 
                1, 1, -1, 
                -1, 1, 1), nrow = 3, ncol = 3, byrow = T)

sign_check(IRF = irf_temp, r_start = 1, r_end = 1, target = 2, K = K, Smat = Smat)


Rcpp::sourceCpp("sign_fast.cpp")

iter = 10
Accept_model = vector("list", length = 10)
SR_EH(Accept_model, K = K, num_slow = 1, C = C, target = 3, Signmat_0 = Signmat_0, Signmat_r = Smat, 
               r_start = 1, r_end = 1, iter = 10, A_hat = A_hat, n_ahead = 5)



