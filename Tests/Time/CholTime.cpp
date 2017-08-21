#include <iostream>
#include <Utils/Utils.hpp>
#include <Regression/SigmaPoints.hpp>


using namespace std;

#define XDIM 200


int main()
{

    unsigned int samples_gaussian = 1000;
    vec mean = zeros<vec>(2);
    mat cov = zeros<mat>(2,2);
    cov(0,0) = 0.3; cov(0,1) = 0.1;
    cov(1,0) = 0.1; cov(1,1) = 0.1;

    mat samples = mat(2, samples_gaussian, fill::zeros);
    const int ncols = cov.n_cols;
    arma_rng::set_seed(0);

    arma::mat Y     = arma::randn(samples_gaussian, ncols);
    samples         = (arma::repmat(mean, 1, samples_gaussian).t() + Y * arma::chol(cov)).t();

    mat save = samples.t();
    save.save("Samples.txt", arma::raw_ascii);

    EllipsoidPoints<2>e_points(samples_gaussian);

    mat eSamples;
    e_points.getSamples(mean, cov, eSamples);

    eSamples.save("SamplesE.txt", arma::raw_ascii);


    return 0;
}
