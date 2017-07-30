#include <iostream>
#include <fstream>

#include <Cost/QuadraticCost.hpp>

using namespace std;

#define XDIM 3
#define UDIM 2


using namespace std;
using namespace arma;



int main()
{

    mat A(XDIM, XDIM, fill::zeros);
    vec a(XDIM, fill::zeros);
    vec in(XDIM, fill::randu);
    double scalar = 0.0;

    A.col(0) = vec(XDIM, fill::ones);
    A.col(1).fill(2.0);

    QuadraticCost<XDIM>test;

    std::cout<<"Value: "        << endl << in <<endl;
    std::cout<<"Evaluation: "   << test.evaluate(in)<<endl;
    std::cout<<"Gradient: "     << endl << test.gradient(in)<<endl;
    std::cout<<"Hessian: "      << endl << test.hessian()<<endl;

    test.quadratize(in, A, a, scalar);

    return 0;

}
