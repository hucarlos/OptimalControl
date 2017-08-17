#include <iostream>
#include <Utils/Utils.hpp>


using namespace std;

#define XDIM 200


int main()
{

    TimeVar  time;
    
    mat X = randu<mat>(XDIM,XDIM);
    mat Y = X.t()*X;
    
    
    time=timeNow();
    mat R2 = chol(Y, "lower");
    double timeex = duration(timeNow() - time);
    
    std::cerr<<std::setprecision(10)<<timeex * 1.0e-6 << std::endl;

    mat R1 = chol(Y);
    
    return 0;
}
