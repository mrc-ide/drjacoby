#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
SEXP loglike(Rcpp::NumericVector params, Rcpp::List data, Rcpp::List misc) {
  
  // extract parameters
  double mu = params["mu"];
  double sigma = params["sigma"];
  
  // unpack data
  std::vector<double> x = Rcpp::as< std::vector<double> >(data["x"]);
  
  // sum log-likelihood over all data
  double ret = 0.0;
  for (unsigned int i = 0; i < x.size(); ++i) {
    ret += R::dnorm(x[i], mu, sigma, true);
  }
  
  // return as SEXP
  return Rcpp::wrap(ret);
}

// [[Rcpp::export]]
SEXP logprior(Rcpp::NumericVector params, Rcpp::List misc) {
  
  // extract parameters
  double sigma = params["sigma"];
  
  // calculate logprior
  double ret = -log(20.0) + R::dlnorm(sigma, 0.0, 1.0, true);
  
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
