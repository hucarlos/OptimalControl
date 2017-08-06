#include <iostream>
#include <Utils/Utils.hpp>
#include <Regression/SigmaPoints.hpp>


using namespace std;

#define XDIM 2

using namespace arma;

typedef typename vec::fixed<XDIM>State;
typedef typename mat::fixed<XDIM, XDIM>StateMat;

int main()
{


    State mean      = ones<vec>(XDIM);
    State initRadius = zeros<vec>(XDIM);

    initRadius(0) = 1;
    initRadius(1) = 1;

    const int samples_factor = 3;
    const int count_points   = std::ceil((XDIM + 3)/8) + 1;
    const int samples4Sigma  = 4*XDIM;

    const int count_samples  = samples_factor * count_points;

    arma::mat points;

    arma::mat allPoints(XDIM, samples4Sigma*count_samples, fill::zeros);


    State radius = initRadius;

    for(int i=0; i<count_samples; i++)
    {


        sigmaPoints<XDIM>(mean, radius, points);

        allPoints.submat(0, i*samples4Sigma, size(XDIM, samples4Sigma)) = points.submat(0,1,size(XDIM, samples4Sigma));

        radius *= 0.25;

    }



    allPoints.save("Points.txt", arma::raw_ascii);

    return 0;
}
