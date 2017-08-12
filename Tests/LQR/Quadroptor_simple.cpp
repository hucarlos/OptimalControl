#include <iostream>
#include <fstream>

#include <Cost/SystemCost.hpp>
#include <System/Quadrotor.hpp>
#include <LQR/iQRLQR.hpp>
#include <LQR/iQRSELQR.hpp>

#include <Utils/SystemParameters.hpp>
#include <Cost/QuadObstacleCost.hpp>

using namespace std;

#define XDIM 12
#define ODIM 3
#define UDIM 4


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

    Quadrotor robot(1.0/20.0);


    const unsigned int ell      = 150;
    const double delta          = 1.0e-4;
    const unsigned int max_iter = 100;

    State xStart          = zeros<mat>(XDIM);
    State xGoal           = zeros<mat>(XDIM);

    // Run iLQR and Extended LQR
    time_t seed = 1372474623; //time_t seed = time(0);
    srand(seed);

    xGoal = zeros<vec>(XDIM);
    int random = rand() % 12;
    xGoal(random / 4)               = ((double) rand() / RAND_MAX) * 4.5 - 2.25; //M_PI; //0; //
    xGoal(((random / 4) + 1) % 3)   = ((random % 4) / 2 == 0 ? 2.25 : -2.25) + ((double) rand() / RAND_MAX) * 0.02 - 0.01;
    xGoal(((random / 4) + 2) % 3)   = ((random % 4) % 2 == 0 ? 2.25 : -2.25) + ((double) rand() / RAND_MAX) * 0.02 - 0.01;

    xStart = -xGoal;
    
    const ControlMat R          = 20  * eye<mat>(UDIM, UDIM);
    const StateMat Q            = 500 * eye<mat>(XDIM, XDIM);


    Control uNominal      = zeros<vec>(UDIM);

    uNominal[0] = uNominal[1] = uNominal[2] = uNominal[3] = robot.getGravity() * robot.getMass()/4;

//    State result = robot.move(xStart, uNominal);


    const double obstacleFactor = 1.0;
    const double scaleFactor    = 10.0;
    const double robotRadius    = 0.3429/2 + 0.1;

    const string mapfile        = "Map1.yaml";

    // For saving results
    std::vector<State>systemPath(ell +1);
    std::vector<Control>nominalControls(ell);
    std::string filename;
    const std::string ext(".txt");

    //============================ Set the obstacles =========================
    vec::fixed<ODIM> bottomLeft, topRight;
    std::vector<Obstacle<ODIM> >obstacles;
    loadMapYAML(mapfile, obstacles, bottomLeft, topRight);

    QuadObstacleCost<XDIM, ODIM>obstacles_cost(robotRadius, obstacles);

    obstacles_cost.setTopRight(topRight);
    obstacles_cost.setBottomLeft(bottomLeft);

    obstacles_cost.setScaleFactor(scaleFactor);
    obstacles_cost.setObstacleFactor(obstacleFactor);

    // =========================== Init all the cost functions ==========================

    std::vector<Control>lNominal(ell);
    std::fill(lNominal.begin(), lNominal.end(), uNominal);

    QuadraticCost<XDIM>init_cost(xStart, Q);
    QuadraticCost<XDIM>final_cost(xGoal, Q);
    QuadraticCost<UDIM>control_cost(uNominal, R);

    SystemCost<XDIM, UDIM, ODIM>system_cost(&control_cost, &obstacles_cost);

    // ========================================= SELQR ALGORITHMS =============================

    SELQR<XDIM, UDIM>selqr(ell, &robot, &init_cost, &system_cost, &final_cost, true);
    t1=timeNow();
    selqr.estimate(xStart, max_iter, delta, lNominal);
    std::cout<<"SELQR TIME(ms) "<<duration(timeNow() - t1)<<std::endl<<endl<<endl;

//    filename = selqr.getName() + ext;
//    selqr.estimatePath(xStart);

//    selqr.getNominalState(systemPath);
//    selqr.getNominalControl(nominalControls);

//    printPathControls<XDIM, UDIM>(systemPath, nominalControls, filename);

    // ========================================= iQRLQR ALGORITHMS =============================

//    iQRLQR<XDIM, UDIM>iqrlqr(ell, &robot, &init_cost, &system_cost, &final_cost, true);
//    iqrlqr.setInitRadius(params.iQRLQR.initRadius);
//    iqrlqr.setEpsilon(params.iQRLQR.epsilon);

//    t1=timeNow();
//    iqrlqr.estimate(xStart, max_iter, delta, lNominal);
//    std::cout<<"iQRLQR TIME(ms) "<<duration(timeNow() - t1)<<std::endl<<endl;

//    filename = iqrlqr.getName() + ext;
//    iqrlqr.estimatePath(xStart);

//    iqrlqr.getNominalState(systemPath);
//    iqrlqr.getNominalControl(nominalControls);

//    printPathControls<XDIM, UDIM>(systemPath, nominalControls, filename);


    // ========================================= iQRSELQR ALGORITHMS =============================

//    iQRSELQR<XDIM, UDIM>qrselqr(ell, &robot, &init_cost, &system_cost, &final_cost, true);
//    qrselqr.setInitRadius(params.iQRSELQR.initRadius);
//    qrselqr.setEpsilon(params.iQRSELQR.epsilon);

//    t1=timeNow();
//    qrselqr.estimate(xStart, max_iter, delta, lNominal);
//    std::cout<<"iQRSELQR TIME(ms) "<<duration(timeNow() - t1)<<std::endl;

//    filename = qrselqr.getName() + ext;
//    std::cout<<"Final cost: "<<qrselqr.estimatePath(xStart)<<endl;

//    qrselqr.getNominalState(systemPath);
//    qrselqr.getNominalControl(nominalControls);

//    printPathControls<XDIM, UDIM>(systemPath, nominalControls, filename);



    return 0;

}
