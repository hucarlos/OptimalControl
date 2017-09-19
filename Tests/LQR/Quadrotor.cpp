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

ExtenedState setDecreceFactors(const double&p, const double&v,
                               const double&r, const double&w,
                               const double&u)
{
    ExtenedState result;
    result.subvec(0, 2)     = p  * ones<vec>(3);
    result.subvec(3, 5)     = v  * ones<vec>(3);
    result.subvec(6, 8)     = r  * ones<vec>(3);
    result.subvec(9, 11)    = w  * ones<vec>(3);
    result.subvec(12, 15)   = u  * ones<vec>(4);

    return result;
}


int main(int argc, char *argv[])
{
    TimeVar t1;

    Quadrotor robot(1.0/25.0);


    const unsigned int ell      = 200;
    const double delta          = 1.0e-4;
    const unsigned int max_iter = 100;

    const double obstacleFactor = 1;
    const double scaleFactor    = 10.0;
    const double robotRadius    = 0.3429/2 + 0.1;

    bool vis = true;

    //============================ Set the obstacles =========================
    const string mapfile        = "map1.yaml";
    vec::fixed<ODIM> bottomLeft, topRight;
    std::vector<Obstacle<ODIM> >obstacles;
    loadMapYAML(mapfile, obstacles, bottomLeft, topRight);

    QuadObstacleCost<XDIM, ODIM>obstacles_cost(robotRadius, obstacles);

    obstacles_cost.setTopRight(topRight);
    obstacles_cost.setBottomLeft(bottomLeft);

    obstacles_cost.setScaleFactor(scaleFactor);
    obstacles_cost.setObstacleFactor(obstacleFactor);


    time_t seed =  /*1503903180;*/time(0);
    srand(seed);

    std::cout<<"Seed: "<<seed<<endl<<endl;

    State xStart          = zeros<mat>(XDIM);
    State xGoal           = zeros<mat>(XDIM);

    int random = rand() % 12;

    xGoal(random / 4)               = ((double) rand() / RAND_MAX) * 4.5 - 2.25; //M_PI; //0; //
    xGoal(((random / 4) + 1) % 3)   = ((random % 4) / 2 == 0 ? 2.25 : -2.25) + ((double) rand() / RAND_MAX) * 0.02 - 0.01;
    xGoal(((random / 4) + 2) % 3)   = ((random % 4) % 2 == 0 ? 2.25 : -2.25) + ((double) rand() / RAND_MAX) * 0.02 - 0.01;

    //        xGoal(0) = 2.2534;
    //        xGoal(1) = 2.2401;
    //        xGoal(2) = -2.1055;

    xStart = -xGoal;

    const ControlMat R  = 15.0  * eye<mat>(UDIM, UDIM);
    const StateMat Q    = 500.0 * eye<mat>(XDIM, XDIM);


    Control uNominal    = zeros<vec>(UDIM);

    uNominal[0] = uNominal[1] = uNominal[2] = uNominal[3] = robot.getGravity() * robot.getMass()/4;

    // For saving results
    std::vector<State>systemPath(ell +1);
    std::vector<Control>nominalControls(ell);
    std::string filename;
    const std::string ext =  (".txt");

    // =========================== Init all the cost functions ==========================

    std::vector<Control>lNominal(ell);
    std::fill(lNominal.begin(), lNominal.end(), uNominal);

    QuadraticCost<XDIM>init_cost(xStart, Q);
    QuadraticCost<XDIM>final_cost(xGoal, Q);
    QuadraticCost<UDIM>control_cost(uNominal, R);

    SystemCost<XDIM, UDIM, ODIM>system_cost(&control_cost, &obstacles_cost);

    try
    {
        ExtenedState initRadius;
        initRadius.subvec(0, 2)     = 1.0e-3  * ones<vec>(3);
        initRadius.subvec(3, 5)     = 1.0e-3 * ones<vec>(3);
        initRadius.subvec(6, 8)     = 1.0e-4  * ones<vec>(3);
        initRadius.subvec(9, 11)    = 1.0e-3  * ones<vec>(3);
        initRadius.subvec(12, 15)   = 1.0e-3  * ones<vec>(4);

        double epsilon      = 1.0e-2;

        // ========================================= SELQR ALGORITHMS =============================

        SELQR<XDIM, UDIM>selqr(ell, &robot, &init_cost, &system_cost, &final_cost, vis);

        t1=timeNow();
        selqr.estimate(xStart, max_iter, delta, lNominal);
        double timeSELQR = duration(timeNow() - t1)/1000.0;

        filename = selqr.getName() + ext;
        selqr.estimatePath(xStart);

        selqr.getNominalState(systemPath);
        selqr.getNominalControl(nominalControls);


        printPathControls<XDIM, UDIM>(systemPath, nominalControls, filename);

        arma::mat out;
        std::string armaName = selqr.getName() + "arma" + ext;
        vector2ArmaMat(systemPath, nominalControls, out);
        out.save(armaName, arma::raw_ascii);


        // ========================================= iQRSELQR ALGORITHMS =============================

        iQRSELQR<XDIM, UDIM>qrselqr(ell, &robot, &init_cost, &system_cost, &final_cost, vis);

        qrselqr.setInitRadius(initRadius);
        qrselqr.setEpsilon(epsilon);

        ExtenedState decrese = setDecreceFactors(0.95, 0.95, 0.95, 0.95, 0.95);
        qrselqr.setDecreceFactors(decrese);
        qrselqr.setMinEig(0.0);
        qrselqr.setFactEig(0.1);

        // Regression parameters
        qrselqr.setSamplingMode(SAMPLING_MODE::ELLIPSOID_S);
        qrselqr.setSamplingFactor(2);
        qrselqr.setParallel(true);

        t1=timeNow();
        qrselqr.estimate(xStart, max_iter, delta, lNominal);
        double timeQRSELQR = duration(timeNow() - t1)/1000.0;

        filename = qrselqr.getName() + ext;
        armaName = qrselqr.getName() + "arma" + ext;

        qrselqr.getNominalState(systemPath);
        qrselqr.getNominalControl(nominalControls);

        vector2ArmaMat(systemPath, nominalControls, out);
        out.save(armaName, arma::raw_ascii);

        printPathControls<XDIM, UDIM>(systemPath, nominalControls, filename);


    }
    catch(std::exception&e)
    {
        std::cerr<<e.what()<<endl;
    }


    return 0;

}
