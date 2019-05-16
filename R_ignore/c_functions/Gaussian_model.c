
#include <math.h>

double loglike(const double *theta, const double *x, int n) {
  
  // extract parameters from theta
  double mu = theta[0];
  double sigma = theta[1];
  
  // sum log-likelihood over all data
  double ret = 0.0;
  for (int i = 0; i < n; ++i) {
    ret += -0.5*log(2*M_PI*sigma*sigma) - (x[i] - mu)*(x[i] - mu)/(2*sigma*sigma);
  }
  return ret;
}

double logprior(const double *theta) {
  
  // extract parameters from theta
  double mu = theta[0];
  double sigma = theta[1];
  
  // calculate logprior
  double ret = -0.5*log(2*M_PI*10*10) - (mu - 0.0)*(mu - 0.0)/(2*10*10);
  ret += -log(sigma) -0.5*log(2*M_PI*1*1) - (log(sigma) - 0)*(log(sigma) - 0)/(2*1*1);
  
  return ret;
}
