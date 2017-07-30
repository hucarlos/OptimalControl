#include <iostream>
#include <fstream>

#include <Cost/ObstacleCost.hpp>

using namespace std;

#define XDIM 3
#define ODIM 2
#define UDIM 2


using namespace std;
using namespace arma;



int main(int argc, char *argv[])
{


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

    arma::vec state(XDIM, fill::zeros);
    state(0) = 0.0;
    state(1) = -25.0;
    state(2) = M_PI;

    std::cout<<"Evaluation: "<<obstacles_cost.evaluate(state)<<endl<<endl;

    mat Q(XDIM, XDIM, fill::zeros);
    vec q(XDIM, fill::zeros);
    double sq;

    obstacles_cost.quadratize(state, Q, q, sq);

    std::cout<<"Q: "<<endl  <<Q <<endl;
    std::cout<<"q: "<<endl  <<q <<endl;

    return 0;

}
