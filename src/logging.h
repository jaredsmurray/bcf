#ifndef GUARD_logging_h
#define GUARD_logging_h
#include <RcppArmadillo.h>


class Logger
{
    private:

    int level;
    int depth;

    public:
    Logger();

    void log(std::string text);
    void setLevel(int inLevel);
    void startContext();
    void stopContext();
    void getVectorHead(Rcpp::NumericVector x, char s[100]);
    void getVectorHead(std::vector<double> x, char s[100]);
    void getVectorHead(double* x, char s[100]);


};

#endif