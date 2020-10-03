#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins("cpp11")]]

using namespace Rcpp;

// matrix exponentials
arma::mat matexp(arma::mat X, int n)
{
  if (n == 0)
  {
    return arma::eye(X.n_cols, X.n_rows);
  }
  else if (n == 1)
  {
    return X;
  }
  else
  {
    return X * matexp(X, n - 1);
  }
}

// IRFs
// [[Rcpp::export]]
List IRF_fast(arma::mat &A_hat, arma::mat &B_hat, int &horizon)
{
  int K = A_hat.n_rows;
  int p = A_hat.n_cols / K;
  List Out(horizon);
  Out[0] = B_hat;
  if (p == 1)
  {
    for (int i = 1; i < horizon; i++)
    {
      Out[i] = matexp(A_hat, i) * B_hat;
    }
    return Out;
  }
  else
  {
    arma::mat Mm(K * p, K * p, arma::fill::zeros);
    Mm.submat(0, 0, (K - 1), (K * p - 1)) = A_hat;
    Mm.submat(K, 0, (K * p - 1), ((p - 1) * K - 1)) = arma::eye(K * (p - 1), K * (p - 1));
    arma::mat Mm1(K * p, K * p, arma::fill::eye);
    for (int i = 0; i < (horizon - 1); i++)
    {
      Mm1 = Mm1 * Mm;
      Out[i + 1] = Mm1.submat(0, 0, (K - 1), (K - 1)) * B_hat;
    }
    return Out;
  }
}

// choose 2 int from N (as index)
arma::imat choose2_fast(int &N)
{
  arma::imat Out(N * (N - 1) / 2, 2);
  int c = 0;
  for (int i = 0; i < (N - 1); i++)
  {
    for (int j = i + 1; j < N; j++)
    {
      Out(c, 0) = i;
      Out(c, 1) = j;
      c++;
    }
  }
  return Out;
}

// Givens rotation
// [[Rcpp::export]]
arma::mat givensQ_fast(arma::vec &thetas, int K)
{
  arma::mat Out = arma::eye(K, K);
  arma::imat Cmat = choose2_fast(K);
  for (int i = 0; i < K * (K - 1) / 2; i++)
  {

    arma::mat temp = arma::eye(K, K);
    temp(Cmat(i, 0), Cmat(i, 0)) = cos(thetas(i));
    temp(Cmat(i, 1), Cmat(i, 1)) = cos(thetas(i));
    temp(Cmat(i, 0), Cmat(i, 1)) = -sin(thetas(i));
    temp(Cmat(i, 1), Cmat(i, 0)) = sin(thetas(i));

    Out = Out * temp;
  }
  return Out.t();
}

// sign check
int sign_check(List &IRF, int &r_start,
               int &r_end, int &target, int &K, arma::mat &Smat)
{
  int h = r_start;
  int Out = 1;

  while (h <= r_end)
  {
    arma::mat irf_h = as<arma::mat>(IRF[h]);
    if (all(all(arma::sign(irf_h.submat(0, target, (K - 1), target)) == Smat.submat(0, target, (K - 1), target))))
    {
      h++;
      Out = 1;
    }
    else
    {
      Out = 0;
      break;
    }
  }

  return Out;
}

// main
// [[Rcpp::export]]
void SR_EH(List & Accept_model, int K, int num_slow, arma::mat & C, int target,
           arma::mat & Signmat_0, arma::mat & Signmat_r, int & r_start,
           int & r_end, int & iter, arma::mat & A_hat, int & n_ahead)
{

  int i = 0;
  //List Accept_model(iter);
  while (i < iter)
  {
    // while loop
    arma::vec thetas = arma::randu<arma::vec>((K - num_slow) * (K - num_slow - 1) / 2);
    thetas = thetas*arma::datum::pi;
    //Rcout << " thetas = " << thetas << std::endl;
    arma::mat R = givensQ_fast(thetas, K - num_slow);
    arma::mat Rotmat = arma::eye(K, K);
    Rotmat.submat(num_slow, num_slow, K - 1, K - 1) = R;
    arma::mat Bmat = arma::inv(Rotmat * C);
  
    for(int m = 0; m < K; m ++)
    {
      if (Bmat(m,m) < 0)
      {
        Bmat.submat(0, m, (K-1), m) = Bmat.submat(0, m, (K-1), m) * -1;
      }
    }
    
    List irf_temp = IRF_fast(A_hat, Bmat, n_ahead);
     
    if (all(all(arma::sign(Bmat) == Signmat_0)) == 0)
    {
      continue;
    }
    else
    {
     
      Rcout << " i = " << i << std::endl;
      
      int h = r_start;
      int check = 1;
      
      while (h <= r_end)
      {
        arma::mat irf_h = as<arma::mat>(irf_temp[h]);
       
        if (all(all(arma::sign(irf_h.submat(0, (target-1), (K - 1), (target-1))) == Signmat_r.submat(0, (target-1), (K - 1), (target-1)))))
        {
          h++;
          check = 1;
        }
        else
        {
          check = 0;
          break;
        }
      }
      
      
      
      if (check)
      {
        Accept_model[i] = irf_temp;
        i++;
      }
      else
      {
        continue;
      }
    }
  }
  
}

