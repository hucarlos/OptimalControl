#include <iostream>
#include <Regression/QuadraticRegression.hpp>
#include <Utils/Utils.hpp>

#include <chrono>
#include <fstream>

using namespace std;
using namespace std;
using namespace arma;

#define XDIM 2
#define UDIM 2

typedef typename vec::fixed<XDIM>State;
typedef typename vec::fixed<XDIM>Control;

typedef typename mat::fixed<XDIM, XDIM>StateMat;
typedef typename mat::fixed<UDIM, UDIM>ControlMat;

typedef typename vec::fixed<XDIM + UDIM>ExtendedState;
typedef typename mat::fixed<XDIM + UDIM, XDIM + UDIM>ExtendedStateMat;

StateMat A      = randn<mat>(XDIM, XDIM);
ControlMat B    = randn<mat>(UDIM, UDIM);

StateMat S      = randn<mat>(XDIM, XDIM);
State s         = randn<vec>(XDIM);

ControlMat R    = randn<mat>(UDIM, UDIM);
StateMat   Q    = randn<mat>(XDIM, XDIM);



/**
 * @brief move
 * @param imput
 * @return
 */
State move(const ExtendedState&imput)
{
    int sub = XDIM + UDIM - 1;

    State x     = imput.subvec(0, XDIM - 1);
    Control u   = imput.subvec(XDIM, sub);

    return A*x + B*u;
}

/**
 * @brief costFunction
 * @param imput
 * @return
 */
double costFunction(const ExtendedState&imput)
{

    int sub = XDIM + UDIM - 1;

    const State     x = imput.subvec(0, XDIM - 1);
    const Control   u = imput.subvec(XDIM, sub);

    const State new_state = move(imput);

    const double instan_cost = 0.5*(dot(u,  R * u)) + 0.5*(dot(x, Q * x));

    const double cost_to_go  = 0.5*(dot(new_state, S * new_state)) + dot(new_state, s);

    return instan_cost + cost_to_go;
}



int main()
{
    ExtendedState xStart    = randn<vec>(XDIM + UDIM);

    //Make symmetric the matrix
    S = symmatu(S);
    Q = symmatu(Q);
    R = symmatu(R);

    std::function<double (const ExtendedState&)>f_cost(costFunction);

    //Create the compose matrix
    ExtendedStateMat compose;

    compose.submat(0,0,size(Q))          = Q + (A.t() * S * A);
    compose.submat(XDIM, XDIM, size(R))  = R + (B.t() * S * B);

    compose.submat(0,XDIM,  size(XDIM, UDIM))     = (A.t() * S * B);
    compose.submat(XDIM, 0, size(UDIM, XDIM))     = (A.t() * S * B).t();

    ExtendedState compose_s = zeros<vec>(XDIM + UDIM);

    int sub = XDIM + UDIM - 1;
    compose_s.subvec(0, XDIM-1)     = A.t() * s;
    compose_s.subvec(XDIM, sub)     = B.t() * s;

    const double alpha                  = 1.0e-2;

    QuadraticRegression<XDIM + UDIM>qr(false);

    // Results matrix
    ExtendedStateMat M;
    ExtendedState m;
    double scalar;

    qr.getRegression(f_cost,xStart, M, m, scalar, alpha, 0.0);

    m -= M*xStart;

    ExtendedStateMat M_diff = (M - compose);
    ExtendedState m_diff    = (m - compose_s);

//    toZero<XDIM + UDIM>(M_diff);
//    toZero<XDIM + UDIM, 1>(m_diff);

    std::cout<<"Diff: "<<endl<<M_diff<<endl;
    std::cout<<"diff: "<<endl<<m_diff<<endl;



    return 0;
}
