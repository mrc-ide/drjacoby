#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
SEXP loglike(Rcpp::NumericVector params, Rcpp::List data, Rcpp::List misc) {
  double mu = params["mu"];
  double gamma = params["gamma"];
  double ret = -gamma*(mu*mu - 1.0)*(mu*mu - 1.0);
  
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
