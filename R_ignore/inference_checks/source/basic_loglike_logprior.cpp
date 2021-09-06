#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
SEXP loglike(Rcpp::NumericVector params, Rcpp::List data, Rcpp::List misc) {
  
  std::vector<double> x = Rcpp::as< std::vector<double> >(data["x"]);
  
  double ret = 0.0;
  for (unsigned int i = 0; i < x.size(); ++i) {
    ret += R::dnorm(x[i], params["mu"], 1.0, true);
  }
  
  return Rcpp::wrap(ret);
}

// [[Rcpp::export]]
SEXP logprior(Rcpp::NumericVector params, Rcpp::List misc) {
  return Rcpp::wrap(0.0);
}

// [[Rcpp::export]]  
SEXP create_xptr(std::string function_name) {  
  typedef SEXP (*funcPtr_likelihood)(Rcpp::NumericVector params, Rcpp::List data, Rcpp::List misc);  
  typedef SEXP (*funcPtr_prior)(Rcpp::NumericVector params, Rcpp::List misc);  
  
  if (function_name == "loglike"){
    return(Rcpp::XPtr<funcPtr_likelihood>(new funcPtr_likelihood(&loglike)));
  } 
  
  if (function_name == "logprior"){
    return(Rcpp::XPtr<funcPtr_prior>(new funcPtr_prior(&logprior)));
  } 

  stop("cpp function %i not found", function_name);
}


