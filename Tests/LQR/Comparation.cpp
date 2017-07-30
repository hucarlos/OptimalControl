#include <iostream>
#include <fstream>

#include <Cost/SystemCost.hpp>
#include <System/DDR.hpp>
#include <LQR/iQRLQR.hpp>
#include <LQR/iQRSELQR.hpp>

#include<Utils/SystemParameters.hpp>

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

    DDR robot(1.0/6.0);

    if(argc<2)
    {
        std::cerr<<"Wrond parameters number"<<endl;
        return -1;
    }

    std::string file(argv[1]);

    ObstaclesParameters<XDIM, UDIM>params;
    ObstacleSystemParameters<XDIM, UDIM>reader;
    reader.loadYAMLParams(file, params);

    const unsigned int ell      = params.horizont;
    const double delta          = params.epsilon;
    const unsigned int max_iter = params.maxIters;

    const State xStart          = params.initState;
    const State xGoal           = params.goalState;

    const ControlMat R          = params.controlFactor * eye<mat>(UDIM, UDIM);
    const StateMat Q            = params.stateFactor   * eye<mat>(XDIM, XDIM);
    const Control uNominal      = params.nominalControl;


    const double obstacleFactor = params.obstacleFactor;
    const double scaleFactor    = params.scaleFactor;
    const double robotRadius    = params.robotRadius;

    const string mapfile        = params.mapName;



    //============================ Set the obstacles =========================
    arma::vec2 bottomLeft, topRight;
    std::vector<Obstacle<ODIM> >obstacles;
    loadMapYAML(mapfile, obstacles, bottomLeft, topRight);

    ObstaclesCost<XDIM, ODIM>obstacles_cost(robotRadius, obstacles);

    obstacles_cost.setTopRight(topRight);
    obstacles_cost.setBottomLeft(bottomLeft);

    obstacles_cost.setScaleFactor(scaleFactor);
    obstacles_cost.setObstacleFactor(obstacleFactor);

    // =========================== Init all the cost functions ==========================

    std::vector<Control>lNominal(ell);
    std::fill(lNominal.begin(), lNominal.end(), zeros<vec>(UDIM));

    QuadraticCost<XDIM>init_cost(xStart, Q);
    QuadraticCost<XDIM>final_cost(xGoal, Q);
    QuadraticCost<UDIM>control_cost(uNominal, R);

    SystemCost<XDIM, UDIM>system_cost(&control_cost, &obstacles_cost);

    // ========================================= SELQR ALGORITHMS =============================

    SELQR<XDIM, UDIM>selqr(ell, &robot, &init_cost, &system_cost, &final_cost, true);
    t1=timeNow();
    selqr.estimate(xStart, max_iter, delta, lNominal);
    std::cout<<"SELQR TIME(ms) "<<duration(timeNow() - t1)<<std::endl<<endl<<endl;

    // ========================================= iQRLQR ALGORITHMS =============================

    iQRLQR<XDIM, UDIM>iqrlqr(ell, &robot, &init_cost, &system_cost, &final_cost, true);
    iqrlqr.setInitRadius(params.iQRLQR.initRadius);
    iqrlqr.setEpsilon(params.iQRLQR.epsilon);

    t1=timeNow();
    iqrlqr.estimate(xStart, max_iter, delta, lNominal);
    std::cout<<"iQRLQR TIME(ms) "<<duration(timeNow() - t1)<<std::endl<<endl;

    // ========================================= iQRSELQR ALGORITHMS =============================

    iQRSELQR<XDIM, UDIM>qrselqr(ell, &robot, &init_cost, &system_cost, &final_cost, true);
    qrselqr.setInitRadius(params.iQRSELQR.initRadius);
    qrselqr.setEpsilon(params.iQRSELQR.epsilon);

    t1=timeNow();
    qrselqr.estimate(xStart, max_iter, delta, lNominal);
    std::cout<<"iQRSELQR TIME(ms) "<<duration(timeNow() - t1)<<std::endl;

    std::vector<State>systemPath(ell +1);
    std::vector<Control>nominalControls(ell);

    std::string name_SELQR("SELQR.txt");
    selqr.estimatePath(xStart);

    selqr.getNominalState(systemPath);
    selqr.getNominalControl(nominalControls);

    printPathControls<XDIM, UDIM>(systemPath, nominalControls, name_SELQR);

    std::string name_iQRLQR("iQRLQR.txt");
    iqrlqr.estimatePath(xStart);

    iqrlqr.getNominalState(systemPath);
    iqrlqr.getNominalControl(nominalControls);

    printPathControls<XDIM, UDIM>(systemPath, nominalControls, name_iQRLQR);

    return 0;

}
