#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
SEXP loglike(Rcpp::NumericVector params, Rcpp::List data, Rcpp::List misc) {
  
  int pi = 0;
  std::vector<double> mu(5);
  for (int i = 0; i < 5; ++i) {
    mu[i] = params[pi++];
  }
  
  int block = misc["block"];
  
  double ret = 0.0;
  if (block == 6) {  // likelihood for the global parameters
    for (int i = 0; i < 5; ++i) {
      ret += R::dnorm(mu[i], 0, 1.0, true);
    }
  } else {  // likelihood for each of the 6 data groups
    NumericVector x = data[block - 1];
    for (int i = 0; i < x.size(); ++i) {
      ret += R::dnorm(x[i], mu[block - 1], 1.0, true);
    }
  }
  
  return Rcpp::wrap(ret);
}

// [[Rcpp::export]]
SEXP logprior(Rcpp::NumericVector params, Rcpp::List misc) {
  return Rcpp::wrap(0.0);
}

// [[Rcpp::export]]  
SEXP create_xptr_loglike(std::string fstr) {  
  typedef SEXP (*funcPtr)(Rcpp::NumericVector params, Rcpp::List data, Rcpp::List misc) ;  
  if (fstr == "loglike"){
    return(Rcpp::XPtr<funcPtr>(new funcPtr(&loglike))) ;
  } else {
    return(Rcpp::XPtr<funcPtr>(R_NilValue)) ; 
  }
}

// [[Rcpp::export]]  
SEXP create_xptr_logprior(std::string fstr) {  
  typedef SEXP (*funcPtr)(Rcpp::NumericVector params, Rcpp::List misc) ;  
  if (fstr == "logprior"){
    return(Rcpp::XPtr<funcPtr>(new funcPtr(&logprior))) ;
  } else {
    return(Rcpp::XPtr<funcPtr>(R_NilValue)) ; 
  }
}  
