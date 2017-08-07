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
    TimeVar tGaussian, tSigma;

    int winnerSigma = 0;

    DDR robot(1.0/6.0);

    // Create distributions for sampling
    std::default_random_engine generator(0);

    unsigned int countExperiments = 100;
    unsigned int experiment = 0;

    if(argc >=2)
    {
        countExperiments = std::stod(argv[1]);
    }

    ofstream out("Cords.txt", std::ofstream::out);
    ofstream win("Wins2.txt", std::ofstream::out);

    State xStart, xGoal;

    do
    {

        std::string file;
        try
        {

            const unsigned int ell      = 100;
            const double delta          = 1.0e-4;
            const unsigned int max_iter = 100;


            std::uniform_real_distribution<double> init_x(-20, 20);
            std::uniform_real_distribution<double> init_y(-24, -40);
            std::uniform_real_distribution<double> init_a(-M_PI, M_PI);

            std::uniform_real_distribution<double> final_x(-20, 20);
            std::uniform_real_distribution<double> final_y(25, 40);
            std::uniform_real_distribution<double> final_a(-M_PI, M_PI);

            xStart(0) = init_x(generator);
            xStart(1) = init_y(generator);
            xStart(2) = init_a(generator);

            xGoal(0) = final_x(generator);
            xGoal(1) = final_y(generator);
            xGoal(2) = final_a(generator);

            const ControlMat R          = 1.0    * eye<mat>(UDIM, UDIM);
            const StateMat Q            = 50.0   * eye<mat>(XDIM, XDIM);
            const Control uNominal      = zeros<vec>(UDIM);


            const double obstacleFactor = 1.0;
            const double scaleFactor    = 1.0;
            const double robotRadius    = 1.675;

            const string mapfile("map1.yaml");

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

            // Same radius and epsilon for all examples
            const double epsilon = 1.0e-1;
            vec::fixed<XDIM + UDIM>radius = ones<vec>(XDIM + UDIM);
            radius(0) = 1;
            radius(1) = 1;
            radius(2) = 0.3;
            radius(3) = 1;
            radius(4) = 1;

            // ========================================= SELQR ALGORITHMS =============================

            iQRSELQR<XDIM, UDIM>gaussianPoints(ell, &robot, &init_cost, &system_cost, &final_cost, false);
            gaussianPoints.setInitRadius(radius);
            gaussianPoints.setEpsilon(epsilon);
            gaussianPoints.setGaussiaSampling(true);
            gaussianPoints.setSamplingFactor(2);

            tGaussian=timeNow();
            gaussianPoints.estimate(xStart, max_iter, delta, lNominal);
            double timeGaussian                = duration(timeNow() - tGaussian);

            // ========================================= iQRSELQR ALGORITHMS =============================

            iQRSELQR<XDIM, UDIM>sigmaPoints(ell, &robot, &init_cost, &system_cost, &final_cost, false);
            sigmaPoints.setInitRadius(radius);
            sigmaPoints.setEpsilon(epsilon);
            sigmaPoints.setGaussiaSampling(false);
            sigmaPoints.setSamplingFactor(1);

            tSigma=timeNow();
            sigmaPoints.estimate(xStart, max_iter, delta, lNominal);
            double timeSigma = duration(timeNow() - tSigma);

            std::cout << experiment <<'\t'
                                    <<"Time (ms): "   << timeGaussian                  <<' '   << timeSigma<<'\t'
                                    <<"Cost: "        << gaussianPoints.getAccum()     <<' '   << sigmaPoints.getAccum()<<'\t'
                                    <<"Iters: "       << gaussianPoints.iterations()   <<' '   << sigmaPoints.iterations()<<endl;

            experiment++;

            if(gaussianPoints.getAccum() > sigmaPoints.getAccum())
            {
                winnerSigma++;
            }

        }
        catch(std::logic_error&e)
        {
            std::cout<<"Example: "<<file<<" not solved: "<<e.what()<<endl;
            out<<xStart(0)<<' '<<xStart(1)<<' '<<(180/M_PI)*xStart(2)<<' '<<xGoal(0)<<' '<<xGoal(1)<<' '<<(180/M_PI)*xGoal(2)
              <<' '<<0<<endl;
            continue;
        }
    }
    while(experiment < countExperiments);

    out.close();
    win.close();

    std::cout<<"Win Sigma%: "<<(double(winnerSigma)/countExperiments)*100<<endl;

    return 0;

}
