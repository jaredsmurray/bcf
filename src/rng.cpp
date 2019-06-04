#include <cmath>
#include "rng.h"

  //standard normal, truncated to be >lo
  double rtnormlo0(double lo) {
    double x;
    if(lo<0) {
      x = R::rnorm(0.0, 1.0);
      while(x<lo) x = R::rnorm(0.0, 1.0);
    } else {
      double a = 0.5*(lo + sqrt(lo*lo + 4.0));
      x = R::rexp(1.0/a) + lo;
      double u = R::runif(0.0, 1.0);
      double diff = (x-a);
      double r = exp(-0.5*diff*diff);
      while(u > r) {
        x = R::rexp(1.0/a) + lo;
        u = R::runif(0.0, 1.0);
        diff = (x-a);
        r = exp(-0.5*diff*diff);
      }
    }
    return x;
  }

  double rtnormlo1(double mean, double lo) {
    return mean + rtnormlo0(lo - mean);
  }

  double rtnormlo(double mean, double sd, double lo) {
    double lostar = (lo-mean)/sd;
    return mean + rtnormlo0(lostar)*sd;
  }

  double rtnormhi1(double mean, double lo) {
    return -rtnormlo1(-mean, -lo);
  }

arma::vec rmvnorm_post(arma::vec &m, arma::mat &Phi) {
  arma::mat R = arma::chol(Phi);
  arma::vec mu = arma::solve(arma::trimatu(R), arma::solve(arma::trimatl(R.t()), m));
  //NumericVector z_ = rnorm(n*mu.size());
  //mat z(z_.begin(), mu.size(), n, false, false);
  arma::vec z(mu.size());
  for(size_t p=0; p<mu.size(); ++p) z(p) = R::rnorm(0.0, 1.0);
  arma::vec res = mu + arma::solve(arma::trimatu(arma::chol(Phi)), z);
  return(res);
}
