#include <RcppArmadillo.h>
#include "logging.h"

using namespace Rcpp;

Logger::Logger(){
    level = 0;
    depth = 0;
}

void Logger::setLevel(int levelIn){
    level = levelIn;
}
void Logger::startContext(){
    depth +=1;
}
void Logger::stopContext(){
    depth +=-1;
}

void Logger::getVectorHead(Rcpp::NumericVector x, char s[100]){
    Rprintf(s,"%f, %f, %f, %f, %f, %f, %f, %f, %f, %f ... ", x[0],x[1],x[2],x[3],x[4],x[5],x[6],x[7],x[8],x[9]);
}

void Logger::getVectorHead(std::vector<double> x, char s[100]){
    Rprintf(s,"%f, %f, %f, %f, %f, %f, %f, %f, %f, %f... ", x[0],x[1],x[2],x[3],x[4],x[5],x[6],x[7],x[8],x[9]);
}

void Logger::getVectorHead(double* x, char s[100]){
    Rprintf(s,"%f, %f, %f, %f, %f, %f, %f, %f, %f, %f... ", x[0],x[1],x[2],x[3],x[4],x[5],x[6],x[7],x[8],x[9]);
}


void Logger::log(std::string text){
    if (level > 0){
        for(int didx=0;didx<depth;didx++){
            Rcout << "--";
        }
        if(depth>0){
            Rcout << " ";
        }

        Rcout << text << "\n";
    }
}
