#include <Rcpp.h>
using namespace Rcpp;

// Log likelihood
// [[Rcpp::export]]
SEXP loglike(Rcpp::NumericVector params, Rcpp::List data, Rcpp::List misc) {
  
  // unpack data
  std::vector<double> x = Rcpp::as< std::vector<double> >(data["x"]);
  
  // unpack params
  double mu = params["mu"];
  double sigma = params["sigma"];
  
  // calculate log-likelihood
  double ret = 0.0;
  for (unsigned int i = 0; i < x.size(); ++i) {
    ret += R::dnorm(x[i], mu, sigma, 1);
  }
  
  // catch underflow
  if (!std::isfinite(ret)) {
    const double OVERFLO_DOUBLE = DBL_MAX/100.0;
    ret = -OVERFLO_DOUBLE;
  }
  
  return Rcpp::wrap(ret);
}

// Log prior
// [[Rcpp::export]]
SEXP logprior(Rcpp::NumericVector params, Rcpp::List misc) {
  
  // unpack params
  double mu = params["mu"];
  double sigma = params["sigma"];
  
  // calculate logprior
  double ret = R::dnorm(mu, 6, 0.1, 1) + R::dnorm(sigma, 1, 0.1, 1);
  
  // catch underflow
  if (!std::isfinite(ret)) {
    const double OVERFLO_DOUBLE = DBL_MAX/100.0;
    ret = -OVERFLO_DOUBLE;
  }
  
  // return as SEXP
  return Rcpp::wrap(ret);
}

// Null log likelihood
// [[Rcpp::export]]
SEXP loglikenull(Rcpp::NumericVector params, Rcpp::List data, Rcpp::List misc) {
  return Rcpp::wrap(0.0);
}

// Null log prior
// [[Rcpp::export]]
SEXP logpriornull(Rcpp::NumericVector params, Rcpp::List misc) {
  return Rcpp::wrap(0.0);
}


// [[Rcpp::export]]  
SEXP create_xptr(std::string function_name) {  
  typedef SEXP (*funcPtr_likelihood)(Rcpp::NumericVector params, Rcpp::List data, Rcpp::List misc);  
  typedef SEXP (*funcPtr_prior)(Rcpp::NumericVector params, Rcpp::List misc);  
  
  if (function_name == "loglike"){
    return(Rcpp::XPtr<funcPtr_likelihood>(new funcPtr_likelihood(&loglike)));
  } 
  if (function_name == "loglikenull"){
    return(Rcpp::XPtr<funcPtr_likelihood>(new funcPtr_likelihood(&loglikenull))) ;
  }
  if (function_name == "logprior"){
    return(Rcpp::XPtr<funcPtr_prior>(new funcPtr_prior(&logprior)));
  } 
  if (function_name == "logpriornull"){
    return(Rcpp::XPtr<funcPtr_prior>(new funcPtr_prior(&logpriornull))) ;
  }
  
  stop("cpp function %i not found", function_name);
}
