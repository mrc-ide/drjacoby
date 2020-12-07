#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
// drjacoby_function_start
SEXP loglike(Rcpp::NumericVector params, Rcpp::List data, Rcpp::List misc) {
  
  std::vector<double> x = Rcpp::as< std::vector<double> >(data["x"]);
  
  double ret = 0.0;
  for (unsigned int i = 0; i < x.size(); ++i) {
    ret += R::dnorm(x[i], params["mu"], 1.0, true);
  }
  
  return Rcpp::wrap(ret);
}
// drjacoby_function_end
