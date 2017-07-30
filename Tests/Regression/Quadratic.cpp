#include <iostream>
#include <Regression/QuadraticRegression.hpp>
#include <Utils/Utils.hpp>


using namespace std;

#define XDIM 3
#define UDIM 3

using namespace std;

typedef typename vec::fixed<XDIM>State;
typedef typename vec::fixed<XDIM>Control;

typedef typename mat::fixed<XDIM, XDIM>StateMat;
typedef typename mat::fixed<UDIM, UDIM>ControlMat;

typedef typename vec::fixed<XDIM + UDIM>ExtendedState;
typedef typename mat::fixed<XDIM + UDIM, XDIM + UDIM>ExtendedStateMat;


ExtendedState s = randn<vec>(XDIM + UDIM);

ControlMat R    = randn<mat>(UDIM, UDIM);
StateMat   Q    = randn<mat>(XDIM, XDIM);

/**
 * @brief costFunction
 * @param imput
 * @return
 */
double costFunction(const ExtendedState&imput)
{

    const State     x = imput.subvec(0, XDIM-1);
    const Control   u = imput.subvec(XDIM, XDIM + UDIM -1);

    const double instan_cost = 0.5*(dot(u, R * u)) + 0.5*(dot(x, Q * x)) + dot(imput, s);

    return instan_cost;
}



int main()
{
    ExtendedState xStart = randn<vec>(XDIM + UDIM);

    //Make symmetric the matrix
    Q = symmatu(Q);
    R = symmatu(R);

    std::function<double (const ExtendedState&)>f_cost(costFunction);


    //Create the compose matrix
    mat::fixed<XDIM + UDIM, XDIM + UDIM>compose = zeros<mat>(XDIM + UDIM, XDIM + UDIM);

    int sub = XDIM + UDIM - 1;
    compose.submat(0,0,XDIM-1, XDIM-1)     = Q;
    compose.submat(XDIM, XDIM, sub, sub)   = R;

    const double alpha                  = 1.0e-2;

    QuadraticRegression<XDIM + UDIM>qr;

    // Results matrix
    ExtendedStateMat M;
    ExtendedState m;
    double scalar;


    qr.getRegression(f_cost, xStart, M, m,   scalar, alpha, 0.0);

    ExtendedStateMat diff =  (compose-M);
    toZero<XDIM + UDIM>(diff);
    std::cout<<"Diff: "<<endl<<diff<<endl;

    m -= M*xStart;
    ExtendedState v_diff = (m-s);
    std::cout<<"diff: "<<endl<<v_diff<<endl;


    return 0;
}
