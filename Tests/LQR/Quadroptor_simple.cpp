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

typedef arma::vec::fixed<XDIM + UDIM>ExtenedState;

typedef arma::mat::fixed<UDIM, XDIM>ControlStateMatrix;

typedef typename QuadraticRegression<XDIM + UDIM>::SAMPLING_MODE SAMPLING_MODE;


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


    xGoal = zeros<vec>(XDIM);

//    // Create distributions for sampling
//    const unsigned int seed = 100;
//    std::default_random_engine generator(seed);
//    std::uniform_real_distribution<double> init_x(-2.26, 2.26);
//    std::uniform_real_distribution<double> init_y(-2.26, 2.26);
//    std::uniform_real_distribution<double> init_z(-2.25, 2.25);

//    xGoal(0) = init_x(generator);
//    xGoal(1) = init_y(generator);
//    xGoal(2) = init_z(generator);

//    time_t seed = 1503359178;
    time_t seed = 1503362689; time(0);
    std::cout<<seed<<endl;
    srand(seed);
    int random = rand() % 12;

    xGoal(random / 4)               = ((double) rand() / RAND_MAX) * 4.5 - 2.25; //M_PI; //0; //
    xGoal(((random / 4) + 1) % 3)   = ((random % 4) / 2 == 0 ? 2.25 : -2.25) + ((double) rand() / RAND_MAX) * 0.02 - 0.01;
    xGoal(((random / 4) + 2) % 3)   = ((random % 4) % 2 == 0 ? 2.25 : -2.25) + ((double) rand() / RAND_MAX) * 0.02 - 0.01;

    xStart = -xGoal;
    
    const ControlMat R          = 20  * eye<mat>(UDIM, UDIM);
    const StateMat Q            = 500 * eye<mat>(XDIM, XDIM);


    Control uNominal      = zeros<vec>(UDIM);

    uNominal[0] = uNominal[1] = uNominal[2] = uNominal[3] = robot.getGravity() * robot.getMass()/4;


    const double obstacleFactor = 1;
    const double scaleFactor    = 10.0;
    const double robotRadius    = 0.3429/2 + 0.1;

    const string mapfile        = "/home/cimat/OptimalControl/Configurations/Quadrotor/map1.yaml";

    // For saving results
    std::vector<State>systemPath(ell +1);
    std::vector<Control>nominalControls(ell);
    std::string filename;
    const std::string ext = std::to_string(seed) + (".txt");

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
    
    ExtenedState initRadius;
    initRadius.subvec(0, 2)     = 1.0e-7  * ones<vec>(3);
    initRadius.subvec(3, 5)     = 1.0e-7 * ones<vec>(3);
    initRadius.subvec(6, 8)     = 1.0e-7  * ones<vec>(3);
    initRadius.subvec(9, 11)    = 1.0e-7  * ones<vec>(3);
    initRadius.subvec(12, 15)   = 1.0e-7  * ones<vec>(4);
    
    double epsilon      = 1.0e-3;
    
   
    // ========================================= SELQR ALGORITHMS =============================

    SELQR<XDIM, UDIM>selqr(ell, &robot, &init_cost, &system_cost, &final_cost, true);

    t1=timeNow();
    selqr.estimate(xStart, max_iter, delta, lNominal);
    std::cout<<"SELQR TIME(s) "<<duration(timeNow() - t1)/1000.0<<std::endl;
    std::cout<<"Final cost: "<<selqr.estimatePath(xStart)<<endl<<endl<<endl;

    filename = selqr.getName() + ext;
    selqr.estimatePath(xStart);

    selqr.getNominalState(systemPath);
    selqr.getNominalControl(nominalControls);

    printPathControls<XDIM, UDIM>(systemPath, nominalControls, filename);


    // ========================================= iQRSELQR ALGORITHMS =============================

    iQRSELQR<XDIM, UDIM>qrselqr(ell, &robot, &init_cost, &system_cost, &final_cost, true);

    
    qrselqr.setInitRadius(initRadius);
    qrselqr.setEpsilon(epsilon);
    qrselqr.setSamplingMode(SAMPLING_MODE::ELLIPSOID_S);
    qrselqr.setSamplingFactor(1);
    qrselqr.setDecreceFactors(0.95);
    qrselqr.setMinEig(0.0);
    qrselqr.setFactEig(0.3);
    qrselqr.setParallel(true);
    
    t1=timeNow();
    qrselqr.estimate(xStart, max_iter, delta, lNominal);
    std::cout<<"iQRSELQR TIME(s) "<<duration(timeNow() - t1)/1000.0<<std::endl;

    filename = qrselqr.getName() + ext;
    std::cout<<"Final cost: "<<qrselqr.estimatePath(xStart)<<endl;

    qrselqr.getNominalState(systemPath);
    qrselqr.getNominalControl(nominalControls);

    printPathControls<XDIM, UDIM>(systemPath, nominalControls, filename);



    return 0;

}
