
#include <Rcpp.h>

using namespace Rcpp;



// [[Rcpp::export]]
NumericVector assignAgePen(NumericVector ageInt, NumericVector AgeDiagIndex, IntegerVector status) {
  int n = ageInt.size();
  NumericVector age_pen(n);
  
  for (int i = 0; i < n; i++) {
    age_pen[i] = status[i] == 0 ? floor(ageInt[i]) : floor(AgeDiagIndex[i]);
  }
  
  return age_pen;
}




// [[Rcpp::export]]
NumericVector calculateLikelihood(NumericVector S0, NumericVector S1, NumericVector RR, IntegerVector status) {
  
  int n = S0.size();
  NumericVector likelihood(n);
  
  for (int i = 0; i < n; i++) {
    if (status[i] == 0) { 
      likelihood[i] = (S1[i] / S0[i]); 
    } else { 
      likelihood[i] = (S1[i] / S0[i]) * pow(RR[i], 1);
    }
  }
  
  return likelihood;
}


// [[Rcpp::export]]
NumericVector calculateStatistics(NumericVector ages) {
  NumericVector result(4);
  result[0] = mean(ages); 
  result[1] = sd(ages);   
  result[2] = min(ages);  
  result[3] = max(ages);  
  return result;
}
