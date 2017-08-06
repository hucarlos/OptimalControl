#include <iostream>
#include <Regression/QuadraticRegression.hpp>
#include <Utils/Utils.hpp>


using namespace std;

#define XDIM 5

using namespace std;

typedef typename vec::fixed<XDIM>State;


typedef typename mat::fixed<XDIM, XDIM>StateMat;



State s = randn<vec>(XDIM);

StateMat   Q    = randn<mat>(XDIM, XDIM);

/**
 * @brief costFunction
 * @param imput
 * @return
 */
double costFunction(const State&imput)
{

    const double instan_cost =  0.5*(dot(imput, Q * imput)) + dot(imput, s);

    return instan_cost;
}

int main()
{

    State xStart = zeros<vec>(XDIM);//randn<vec>(XDIM);

    //Make symmetric the matrix
    Q = symmatu(Q);

    std::function<double (const State&)>f_cost(costFunction);


    //Create the compose matrix
    mat::fixed<XDIM, XDIM>compose = zeros<mat>(XDIM, XDIM);

    compose.submat(0,0,XDIM-1, XDIM-1)     = Q;

    const double alpha                  = 1.0e-2;

    QuadraticRegression<XDIM>qr(false, 4);

    // Results matrix
    StateMat M;
    State m;
    double scalar;


    qr.getRegression(f_cost, xStart, M, m,   scalar, alpha, 0.0);

    StateMat diff =  (compose-M);
//    toZero<XDIM, XDIM>(diff);
    std::cout<<"Diff: "<<endl<<diff<<endl;

    m -= M*xStart;
    State v_diff = (m-s);
//    toZero<XDIM, 1>(v_diff);
    std::cout<<"diff: "<<endl<<v_diff<<endl;


    return 0;
}
