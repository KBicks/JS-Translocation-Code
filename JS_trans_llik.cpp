#include <Rcpp.h>
using namespace Rcpp;


// [[Rcpp::export]]
double calc_llik_T(  int K,
                     double n,
                     double N,
                     IntegerVector en, 
                     IntegerVector ex, 
                     NumericVector p, 
                     NumericVector phi, 
                     IntegerMatrix matrix_ch) {  
  
  
  int nrow = matrix_ch.nrow();
  
  NumericVector prob_t(nrow);
  
  for (int i = 0; i < nrow+1; i++) {  
    
    for(int d = ex[i]-1; d < K; d++){
      
      double prod_phi = 1;
      
      if(d != 0){
        for(int j = en[i]-1; j < d; j++){
          prod_phi *= phi[j];
        }
      }
      
      double prod_p = 1;
      
      for(int j = en[i]; j < d+1; j++){
        
        prod_p *=  pow(p[j], matrix_ch(i,j)) * pow(1-p[j], 1-matrix_ch(i,j)); 
      }
      
      
      
      prob_t[i] += (prod_phi * (1-phi[d]) * prod_p);
      
    }
    
    
  }
  
  double L = sum(log(prob_t));
  L *= -1;
  return L;
  
}

// ------------------------------------------------------------------------------------------------------

// [[Rcpp::export]]
double calc_llik_W(int tau,
                int K,
                double n,
                double N,
                IntegerVector en, 
                IntegerVector ex, 
                NumericVector p, 
                NumericVector phi, 
                NumericVector beta, 
                IntegerMatrix matrix_ch) {  
  
  
  int nrow = matrix_ch.nrow();
  
  NumericVector prob_t(nrow);
  
  for (int i = 0; i < nrow+1; i++) {  
    
    for(int b = tau-1; b < en[i]; b++){
      
      for(int d = ex[i]-1; d < K; d++){
        
        double prod_phi = 1;
        
        if(b != d){
          for(int j = b; j < d; j++){
            prod_phi *= phi[j];
          }
          
        }
        
        double prod_p = 1;
        
        for(int j = b; j < d+1; j++){
          
          prod_p *=  pow(p[j], matrix_ch(i,j)) * pow(1-p[j], 1-matrix_ch(i,j)); 
        }
        
        
        
        prob_t[i] += (beta[b] * prod_phi * (1-phi[d]) * prod_p);
        
      }
      
    }
  }
  
  
  
  double prob_0 = 0;
  
  for(int b = tau-1; b < K; b++){
    for(int d = b; d < K; d++){
      
      double prod_phi_un = 1;
      
      if(b != d){
        for(int j = b; j < d; j++){
          prod_phi_un *= phi[j];
        }
      }
      
      double  prod_p_un = 1;
      
      for(int j = b; j < d+1; j++){
        prod_p_un *= (1-p[j]);
      } 
      
      prob_0 +=  (beta[b] * prod_phi_un * (1-phi[d]) * prod_p_un);
      
    }
  }
  
  
  double L = lgamma(N+1) - lgamma(n+1) + sum(log(prob_t)) + n*log(prob_0); 
  L *= -1;
  return L;
  
}

// ------------------------------------------------------------------------------------------------------------------