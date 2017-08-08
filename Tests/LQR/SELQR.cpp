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
    TimeVar tSELQR, tQRSELQR;

    int winner = 0;

    DDR robot(1.0/6.0);

    // Create distributions for sampling
    std::default_random_engine generator(26);

    ofstream out("Cords.txt", std::ofstream::out);
    ofstream win("WinsrLQR.txt", std::ofstream::out);


    unsigned int countExperiments = 100;
    unsigned int experiment = 0;
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
            std::uniform_real_distribution<double> init_a(0, M_PI_2);

            std::uniform_real_distribution<double> final_x(-20, 20);
            std::uniform_real_distribution<double> final_y(25, 40);
            std::uniform_real_distribution<double> final_a(0, M_PI_2);

            State xStart          = zeros<vec>(XDIM);
            xStart(0) = init_x(generator);
            xStart(1) = init_y(generator);
            xStart(2) = init_a(generator);

            State xGoal           = zeros<vec>(XDIM);
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
            const double epsilon = 1.0e-2;
            vec::fixed<XDIM + UDIM>radius = ones<vec>(XDIM + UDIM);
            radius(0) = 1;
            radius(1) = 1;
            radius(2) = 0.1;
            radius(3) = 1;
            radius(4) = 1;

            // ========================================= SELQR ALGORITHMS =============================

            SELQR<XDIM, UDIM>selqr(ell, &robot, &init_cost, &system_cost, &final_cost, false);
            tSELQR=timeNow();
            selqr.estimate(xStart, max_iter, delta, lNominal);
            double timeSELQR = duration(timeNow() - tSELQR);

            // ========================================= iQRSELQR ALGORITHMS =============================

            iQRLQR<XDIM, UDIM>iqrlqr(ell, &robot, &init_cost, &system_cost, &final_cost, false);
            iqrlqr.setInitRadius(radius);
            iqrlqr.setEpsilon(epsilon);
            iqrlqr.setGaussiaSampling(false);
            iqrlqr.setSamplingFactor(2);

            tQRSELQR=timeNow();
            iqrlqr.estimate(xStart, max_iter, delta, lNominal);
            double timeQRSELQR = duration(timeNow() - tQRSELQR);

            if(selqr.getAccum() > iqrlqr.getAccum())
            {
                winner ++;
                win << timeSELQR            <<'\t'   << timeQRSELQR<<'\t'
                    << selqr.getAccum()     <<'\t'   << iqrlqr.getAccum()<<'\t'
                    << selqr.iterations()   <<'\t'   << iqrlqr.iterations()<<'\t'
                    << 0.0                  <<'\t'   << iqrlqr.getEpsilon() << endl;
            }
            else
            {
                out<<xStart(0)<<' '<<xStart(1)<<' '<<(180/M_PI)*xStart(2)<<' '<<xGoal(0)<<' '<<xGoal(1)<<' '<<(180/M_PI)*xGoal(2)<<endl;
            }


            std::cout << experiment <<'\t'
                                    <<"Time (ms): "   << timeSELQR            <<' '   << timeQRSELQR<<'\t'
                                    <<"Cost: "        << selqr.getAccum()     <<' '   << iqrlqr.getAccum()<<'\t'
                                    <<"Iters: "       << selqr.iterations()   <<' '   << iqrlqr.iterations()<<endl;

            experiment++;

        }
        catch(std::logic_error&e)
        {
            std::cout<<"Example: "<<file<<" not solved: "<<e.what()<<endl;
                       continue;
        }
    }
    while(experiment < countExperiments);

    std::cout<<"Win %: "<<(double(winner)/countExperiments)*100<<endl;

    win.close();
    out.close();

    return 0;

}
