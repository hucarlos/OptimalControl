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

typedef arma::vec::fixed<XDIM + UDIM>ExtendedState;

typedef arma::mat::fixed<UDIM, XDIM>ControlStateMatrix;

typedef typename QuadraticRegression<XDIM + UDIM>::SAMPLING_MODE SAMPLING_MODE;


int main(int argc, char *argv[])
{
    TimeVar t1;

    DDR robot(0.1);

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

    // For saving reslts
    std::vector<State>systemPath(ell +1);
    std::vector<Control>nominalControls(ell);
    std::string filename;
    const int index = 8;
    const std::string ext = "Epsilon" + to_string(index) + ".txt";

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

    filename = selqr.getName() + ext;
    selqr.estimatePath(xStart);

    selqr.getNominalState(systemPath);
    selqr.getNominalControl(nominalControls);

    printPathControls<XDIM, UDIM>(systemPath, nominalControls, filename);

    // ========================================= iQRLQR ALGORITHMS =============================

//    iQRLQR<XDIM, UDIM>iqrlqr(ell, &robot, &init_cost, &system_cost, &final_cost, true);
//    iqrlqr.setInitRadius(params.iQRLQR.initRadius);
//    iqrlqr.setEpsilon(params.iQRLQR.epsilon);

//    vec::fixed<XDIM + UDIM>decres1 = 0.5 * ones<vec>(XDIM + UDIM);
//    decres1(2) = 0.1;

//    iqrlqr.setSamplingMode(SAMPLING_MODE::ELLIPSOID_S);
//    iqrlqr.setSamplingFactor(1);
//    iqrlqr.setDecreceFactors(decres1);
//    iqrlqr.setMinEig(0.0);
//    iqrlqr.setFactEig(0.1);

//    t1=timeNow();
//    iqrlqr.estimate(xStart, max_iter, delta, lNominal);
//    std::cout<<"iQRLQR TIME(ms) "<<duration(timeNow() - t1)<<std::endl<<endl;

//    filename = iqrlqr.getName() + ext;
//    iqrlqr.estimatePath(xStart);

//    iqrlqr.getNominalState(systemPath);
//    iqrlqr.getNominalControl(nominalControls);

//    printPathControls<XDIM, UDIM>(systemPath, nominalControls, filename);


    // ========================================= iQRSELQR ALGORITHMS =============================

    iQRSELQR<XDIM, UDIM>qrselqr(ell, &robot, &init_cost, &system_cost, &final_cost, true);
    qrselqr.setInitRadius(params.iQRSELQR.initRadius);
    qrselqr.setEpsilon(params.iQRSELQR.epsilon);

    qrselqr.setSamplingMode(SAMPLING_MODE::ELLIPSOID_S);
    qrselqr.setSamplingFactor(3);

    vec::fixed<XDIM + UDIM>decres2 = 0.75 * ones<vec>(XDIM + UDIM);
    decres2(0) = 0.6;
    decres2(2) = 0.1;

    qrselqr.setDecreceFactors(decres2);
    qrselqr.setMinEig(0.0);
    qrselqr.setFactEig(0.1);
    qrselqr.setParallel(false);

    t1=timeNow();
    qrselqr.estimate(xStart, max_iter, delta, lNominal);
    std::cout<<"iQRSELQR TIME(ms) "<<duration(timeNow() - t1)<<std::endl;

    filename = qrselqr.getName() + ext;
    std::cout<<"Final cost: "<<qrselqr.estimatePath(xStart)<<endl;

    qrselqr.getNominalState(systemPath);
    qrselqr.getNominalControl(nominalControls);

    printPathControls<XDIM, UDIM>(systemPath, nominalControls, filename);

    std::vector< arma::vec::fixed<XDIM + UDIM> >vec_radius;
    qrselqr.getVecRadius(vec_radius);

    string r_name = "Radius" + ext;
    ofstream outradius(r_name);
    for(unsigned int i=0; i<vec_radius.size(); i++)
    {
        ExtendedState rad = vec_radius.at(i);
        for(int j=0; j<rad.n_elem; j++)
        {
            outradius<<rad(j)<<' ';
        }
        outradius<<endl;
    }

    return 0;

}
