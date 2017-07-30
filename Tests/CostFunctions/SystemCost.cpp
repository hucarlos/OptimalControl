#include <iostream>
#include <fstream>

#include <Cost/SystemCost.hpp>

using namespace std;

#define XDIM 3
#define ODIM 2
#define UDIM 2


using namespace std;
using namespace arma;

typedef arma::vec::fixed<XDIM>State;
typedef arma::vec::fixed<UDIM>Control;

typedef arma::mat::fixed<XDIM, XDIM>StateMat;
typedef arma::mat::fixed<UDIM, UDIM>ControlMat;

typedef arma::mat::fixed<XDIM + UDIM, 1>ExtenedState;

typedef arma::mat::fixed<UDIM, XDIM>ControlStateMatrix;


int main(int argc, char *argv[])
{
    arma::vec state(XDIM, fill::zeros);
    state(0) = 0.0;
    state(1) = -25.0;
    state(2) = M_PI;

    std::string mapfile(argv[1]);
    const double obstacleFactor = 1.0;
    const double scaleFactor    = 1.0;
    const double robotRadius    = 3.35/2.0;

    arma::vec2 bottomLeft, topRight;
    std::vector<Obstacle<ODIM> >obstacles;


    loadMapYAML(mapfile, obstacles, bottomLeft, topRight);

    //============================ Set the obstacles =========================
    ObstaclesCost<XDIM, ODIM>obstacles_cost(robotRadius, obstacles);

    obstacles_cost.setTopRight(topRight);
    obstacles_cost.setBottomLeft(bottomLeft);

    obstacles_cost.setScaleFactor(scaleFactor);
    obstacles_cost.setObstacleFactor(obstacleFactor);

    // =========================== Init all the cost functions ==========================

    const ControlMat R(eye<mat>(UDIM, UDIM));
    const StateMat Q(zeros<mat>(XDIM, XDIM));

    State xStart(zeros<mat>(XDIM));
    xStart(0) = 0.0;
    xStart(1) = -25.0;
    xStart(2) = M_PI;

    State xGoal(zeros<vec>(XDIM));

    const Control uNominal(zeros<mat>(UDIM));

    QuadraticCost<XDIM>init_cost(xStart, Q);
    QuadraticCost<XDIM>final_cost(xGoal, Q);
    QuadraticCost<UDIM>control_cost(uNominal, R);
    SystemCost<XDIM, UDIM>system_cost(&control_cost, &obstacles_cost);

    std::cout<<"Evaluation: "<<system_cost.evaluate(xStart, uNominal)<<endl<<endl;

    StateMat H(zeros<mat>(XDIM, XDIM));
    State g(zeros<vec>(XDIM));

    ControlMat Hc(zeros<mat>(UDIM, UDIM));
    Control gc(zeros<vec>(UDIM));

    ControlStateMatrix Pt(zeros<mat>(UDIM, XDIM));

    system_cost.quadratize(xStart, uNominal, H, Hc, Pt, g, gc);

    std::cout<<"Hessian State: "    <<endl  << H    <<endl;
    std::cout<<"Gradient State: "   <<endl  << g    <<endl;

    std::cout<<"Hessian Control: "    <<endl  << Hc    <<endl;
    std::cout<<"Gradient Control: "   <<endl  << gc    <<endl;


    return 0;

}
