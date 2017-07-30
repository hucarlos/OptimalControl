#include <iostream>
#include <fstream>

#include <Cost/SystemCost.hpp>
#include <System/DDR.hpp>
#include <LQR/iQRLQR.hpp>
#include <LQR/iQRSELQR.hpp>

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
    TimeVar t1;

    const unsigned int ell = 100;

    DDR robot(1.0/6.0);

    if(argc<2)
    {
        std::cerr<<"Wrond parameters number"<<endl;
        return -1;
    }

    std::string mapfile(argv[1]);
    const double obstacleFactor = 2.0;
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

    const ControlMat R = 1.0  * eye<mat>(UDIM, UDIM);
    const StateMat   Q = 50.0 * eye<mat>(XDIM, XDIM);

    State xStart;
    xStart(0) = 13;
    xStart(1) = 23.0;
    xStart(2) = 15*(M_PI/180.0);

    State xGoal;
    xGoal(0) = -14;
    xGoal(1) = -7;
    xGoal(2) = 45*(M_PI/180.0);;

    const Control uNominal = zeros<vec>(UDIM);

    std::vector<Control>lNominal(ell);
    std::fill(lNominal.begin(), lNominal.end(), zeros<vec>(UDIM));

    QuadraticCost<XDIM>init_cost(xStart, Q);
    QuadraticCost<XDIM>final_cost(xGoal, Q);
    QuadraticCost<UDIM>control_cost(uNominal, R);

    // Radius search
    ExtenedState radius = ones<vec>(XDIM+UDIM);
    const double epsilon = std::stod(argv[2]);

    SystemCost<XDIM, UDIM>system_cost(&control_cost, &obstacles_cost);

    iQRLQR<XDIM, UDIM>iqrselqr(ell, &robot, &init_cost, &system_cost, &final_cost, true);
    iqrselqr.setInitRadius(radius);
    iqrselqr.setEpsilon(epsilon);


    iQRSELQR<XDIM, UDIM>qrselqr(ell, &robot, &init_cost, &system_cost, &final_cost, true);
    qrselqr.setInitRadius(radius);
    qrselqr.setEpsilon(epsilon);

    t1=timeNow();
    iqrselqr.estimate(xStart, 100, 1.0e-3, lNominal);
    std::cout<<"iQRLQR TIME(ms) "<<duration(timeNow() - t1)<<std::endl<<endl;

    t1=timeNow();
    qrselqr.estimate(xStart, 100, 1.0e-3, lNominal);
    std::cout<<"iQRSELQR TIME(ms) "<<duration(timeNow() - t1)<<std::endl;

    std::vector<State>systemPath(ell +1);
    std::vector<Control>nominalControls(ell);

    std::string name_regre("SELQR.txt");
    std::string filename2 =  name_regre;
    qrselqr.estimatePath(xStart);

    qrselqr.getNominalState(systemPath);
    qrselqr.getNominalControl(nominalControls);


    printPathControls<XDIM, UDIM>(systemPath, nominalControls, filename2);

    return 0;

}
