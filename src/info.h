#ifndef GUARD_info_h
#define GUARD_info_h

#include<vector>

//data
class dinfo {
public:
   dinfo() {p=0;n=0;x=0;y=0;}
   size_t p;  //number of vars
   size_t n;  //number of observations
   double *x; // jth var of ith obs is *(x + p*i+j)
   double *y; // ith y is *(y+i) or y[i]
};

//prior and mcmc
class pinfo
{
public:
   pinfo() {pbd=1.0;pb=.5;alpha=.95;beta=.5;tau=1.0;sigma=1.0;minperbot=5;taccept=0;tavgd=0;tmaxd=0;}
//mcmc info
   double pbd; //prob of birth/death
   double pb;  //prob of birth
//prior info
   double alpha;
   double beta;
   double tau;
//sigma
   double sigma;
//tree min obs
  size_t minperbot;
//mcmc stuffs
  std::vector< std::vector<double> > cvpm;  // Change of variable proposal probability matrix
  unsigned int taccept; //acceptance count of tree proposals
  unsigned int tproposal; //number of tree proposals
  unsigned int tavgd; //average tree depth
  unsigned int tmaxd; //maximum tree depth
};

//sufficient statistics for 1 node
class sinfo
{
public:
   sinfo() {n0=0.0;n=0;sy=0.0;sy2=0.0;}
   double n0; //unweighted sample counts
   double n;
   double sy;
   double sy2;
};

#endif
