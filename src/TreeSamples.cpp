//#include <Rcpp.h>
#include <RcppArmadillo.h>

#include <iostream>
#include <fstream>
#include <vector>
#include <ctime>

#include "rng.h"
#include "tree.h"
#include "info.h"
#include "funs.h"
#include "bd.h"

using namespace Rcpp;


class TreeSamples {
  public:
  bool init;
  size_t m,p,ndraws;
	xinfo xi;
	std::vector<std::vector<tree> > t;


  void load(CharacterVector treef_name_) {
    Rcout << "Loading...\n";
    std::string treef_name = as<std::string>(treef_name_); 
    std::ifstream treef(treef_name.c_str());
	  treef >> xi; //load the cutpoints
    //prxi(xi);
	  treef >> m;  //number of trees
    Rcout << "ntrees " << m << endl;
	  treef >> p;  //dimension of x's
    Rcout << "p " << p << endl;
	  treef >> ndraws; //number of draws from the posterior that were saved.
    Rcout << "ndraws " << ndraws << endl;
    
    t.resize(ndraws,std::vector<tree>(m));
    for(size_t i=0;i<ndraws;i++) {
		  for(size_t j=0;j<m;j++) {
			  treef >> t[i][j];
		  }
    }
    Rcout << "done" <<endl;
    init = true;

    treef.close();
  }
  
  NumericMatrix predict(NumericMatrix x_) {
    size_t n = x_.ncol();
    NumericMatrix ypred(ndraws, n);
    if(init) {
      std::vector<double> x;
      for(NumericMatrix::iterator it=x_.begin(); it!=x_.end(); ++it) x.push_back(*it);
      
      dinfo di;
      di.n=n; di.p=p; di.x = &x[0]; di.y=0;

      for(size_t k=0; k<n; ++k) {
        for(size_t i=0;i<ndraws;i++) {
    		  ypred(i,k) += fit_i(k, t[i], xi, di);
        }
      }
    } else {
      Rcout << "Uninitialized" <<'\n';
    }
    return ypred;
  }
  
  //predictions for multiplicative trees (precision)
  NumericMatrix predict_prec(NumericMatrix x_) {
    size_t n = x_.ncol();
    NumericMatrix ypred(ndraws, n);
    ypred.fill(1.0);
    if(init) {
      std::vector<double> x;
      for(NumericMatrix::iterator it=x_.begin(); it!=x_.end(); ++it) x.push_back(*it);
      
      dinfo di;
      di.n=n; di.p=p; di.x = &x[0]; di.y=0;

      for(size_t k=0; k<n; ++k) {
        for(size_t i=0;i<ndraws;i++) {
      	  ypred(i,k) *= fit_i_mult(k, t[i], xi, di);
        }
      }
    } else {
      Rcout << "Uninitialized" <<'\n';
    }
    return ypred;
  }
  
  //predictions from the ith mcmc iterate
  NumericMatrix predict_i(NumericMatrix x_, size_t i) {
    size_t n = x_.ncol();
    NumericMatrix ypred(1, n);
    if(init) {
      std::vector<double> x;
      for(NumericMatrix::iterator it=x_.begin(); it!=x_.end(); ++it) x.push_back(*it);
      
      dinfo di;
      di.n=n; di.p=p; di.x = &x[0]; di.y=0;

      for(size_t k=0; k<n; ++k) {
        //for(size_t i=0;i<ndraws;i++) {
      	  ypred(0,k) += fit_i(k, t[i], xi, di);
        //}
      }
    } else {
      Rcout << "Uninitialized" <<'\n';
    }
    return ypred;
  }
  
  //predictions from the ith mcmc iterate
  NumericMatrix predict_prec_i(NumericMatrix x_, size_t i) {
    size_t n = x_.ncol();
    NumericMatrix ypred(1, n); ypred.fill(1.0);
    if(init) {
      std::vector<double> x;
      for(NumericMatrix::iterator it=x_.begin(); it!=x_.end(); ++it) x.push_back(*it);
      
      dinfo di;
      di.n=n; di.p=p; di.x = &x[0]; di.y=0;

      for(size_t k=0; k<n; ++k) {
        //for(size_t i=0;i<ndraws;i++) {
          ypred(0,k) *= fit_i_mult(k, t[i], xi, di);
        //}
      }
    } else {
      Rcout << "Uninitialized" <<'\n';
    }
    return ypred;
  }
  
  
  TreeSamples() {}
  
};

RCPP_MODULE(TreeSamples) {
  class_<TreeSamples>( "TreeSamples" )
  .constructor()
  .method( "load", &TreeSamples::load )
  .method( "predict", &TreeSamples::predict  )
  .method( "predict_prec", &TreeSamples::predict_prec  )
  .method( "predict_i", &TreeSamples::predict_i  )
  .method( "predict_prec_i", &TreeSamples::predict_prec_i  )
  ;
}
