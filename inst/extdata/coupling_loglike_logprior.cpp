#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
SEXP loglike(Rcpp::NumericVector params, Rcpp::List data, Rcpp::List misc) {
  
  // unpack data
  std::vector<double> x = Rcpp::as< std::vector<double> >(data["x"]);
  
  // unpack parameters
  double alpha = params["alpha"];
  double beta = params["beta"];
  double epsilon = params["epsilon"];
  
  // sum log-likelihood over all data
  double mean = alpha*alpha*beta + epsilon;
  double ret = 0.0;
  for (unsigned int i = 0; i < x.size(); ++i) {
    ret += R::dnorm(x[i], mean, 1.0, true);
  }
  
  // return as SEXP
  return Rcpp::wrap(ret);
}

// [[Rcpp::export]]
SEXP logprior(Rcpp::NumericVector params, Rcpp::List misc) {
  
  // unpack parameters
  double epsilon = params["epsilon"];
  
  // calculate logprior
  double ret = -log(20.0) - log(10.0) + R::dnorm(epsilon, 0.0, 1.0, true);
  
  // return as SEXP
  return Rcpp::wrap(ret);
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
