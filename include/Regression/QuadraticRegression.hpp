#ifndef QUADRATICREGRESSION_HPP
#define QUADRATICREGRESSION_HPP

#include <iostream>
#include <armadillo>

//#include <Regression/SigmaPoints.hpp>

using namespace std;
using namespace arma;

template<uword xDim>
class QuadraticRegression
{
    public:
        typedef arma::vec::fixed<xDim>Inputs;
        typedef arma::mat::fixed<xDim, xDim>InputsMat;

        QuadraticRegression(): gaussian_sampling(true), samples_factor(2), seed(0)
        {

        }

        QuadraticRegression(const bool gauss, const uint64_t&s): gaussian_sampling(gauss), samples_factor(2), seed(s)
        {

        }

        QuadraticRegression(const bool gauss, const uint64_t&s, const unsigned int sf):
            gaussian_sampling(gauss), samples_factor(sf), seed(s)
        {

        }

    public:

        /**
         * @brief getRegression
         * @param cost
         * @param mean
         * @param A
         * @param a
         * @param scalar
         * @param covar
         * @param lambda
         * @return
         */
        double getRegression(std::function<double (const Inputs&)>cost,
                             const Inputs&mean,
                             InputsMat&A,
                             Inputs&a,
                             double&scalar,
                             const Inputs&radius,
                             const double&lambda=0.0)
        {

            const unsigned int unknow           = ((xDim * (xDim+1))/2) + xDim;
            const unsigned int samples_gaussian = samples_factor * unknow;

            arma::mat samples(xDim, samples_gaussian, fill::zeros);

            InputsMat cov   = zeros<mat>(xDim, xDim);
            cov.diag()      = radius;

            if(gaussian_sampling)
            {

                const int ncols = cov.n_cols;
                arma_rng::set_seed(seed);
                arma::mat Y     = arma::randn(samples_gaussian, ncols);
                samples         = (arma::repmat(mean, 1, samples_gaussian).t() + Y * arma::chol(cov)).t();
            }
            else
            {
//                MatrixXd points_i;
//                const unsigned int count_sigma_point = (2*xDim + 1);
//                const unsigned int samples_sigma = samples_factor * count_sigma_point;

//                samples = MatrixXd::Zero(xDim, samples_sigma);

//                const double startAlpha = 1.0e-2;

//                for(unsigned int si=0; si<samples_factor; si++)
//                {

//                    double alpha = startAlpha * 0.1;
//                    sigmaPoints<xDim>(mean, cov, points_i, alpha);

//                    // Block of size (p,q), starting at (i,j)

//                    int begin_col = si * count_sigma_point;

//                    samples.block(0,begin_col, xDim, count_sigma_point) = points_i;
//                }
            }

            double relative_error = this->solve(cost, mean, samples, A, a, scalar, lambda);
            return relative_error;
        }


                /**
                 * @brief getRegression
                 * @param cost
                 * @param state
                 * @param A
                 * @param a
                 * @param alpha
                 * @param lambda
                 * @param print
                 * @return
                 */
                double getRegression(std::function<double (const Inputs&)>cost,
                                     const Inputs&mean,
                                     InputsMat&A,
                                     Inputs&a,
                                     double&scalar,
                                     const double&alpha=0.01,
                                     const double&lambda=0.0)
                {
                    Inputs radius = alpha * ones<vec>(xDim);
                    return this->getRegression(cost, mean, A, a, scalar, radius, lambda);
                }


        /**
         * @brief solve
         * @param cost
         * @param mean
         * @param samples
         * @param A
         * @param a
         * @param scalar
         * @param lambda
         * @return
         */
        double solve(std::function<double (const Inputs&)>cost,
                         const Inputs&mean,
                         const mat&samples,
                         InputsMat&A,
                         Inputs&a,
                         double&scalar,
                         const double&lambda=0.0)
        {


            const int mat_unknow    = (xDim*(xDim + 1)/2);
            const int unknow        =  mat_unknow + xDim;

            if(unknow <= 0)
            {
                throw std::runtime_error("Bad dimensions for unknow on solve");
            }

            const int samples_count = samples.n_cols;


           mat X(samples_count, unknow, fill::zeros);
           vec y(samples_count, fill::zeros);


            // Cicle for every sample
            const double nominal_value = cost(mean);

            for(int i=0; i<samples_count; i++)
            {
                const Inputs sample = samples.col(i);
                const double value  = cost(sample);
                const Inputs diff   = sample - mean;

                // Cicle for all the matrix values
                int p = 0;
                for(unsigned int r=0; r<xDim; r++)
                {
                    for(unsigned int c=r; c<xDim; c++, p++)
                    {
                        X(i,p)  = diff(r)*diff(c);

                        X(i,p) *= 0.5;

                        if(c>r)
                        {
                            X(i,p) *= 2.0;
                        }
                    }
                }

                for(int r=0; r<xDim; r++, p++)
                {
                    X(i,p) = diff(r);
                }

                y(i) = (value - nominal_value);

            }

            // Obtain the solution
            const mat::fixed<unknow, unknow>I = eye<mat>(unknow, unknow);

            mat::fixed<unknow, unknow> L    = X.t() * X;
            vec::fixed<unknow>b             = X.t() * y;

            mat::fixed<unknow, unknow> M    = (L + lambda*I);
            vec::fixed<unknow> sol;

            sol = arma::solve(M, b, solve_opts::equilibrate);


            // Put results on A
            for(unsigned int i=0, p=0; i<xDim; i++)
            {
                for(unsigned int j=i; j<xDim; j++, p++)
                {
                    A(i,j) = A(j,i) = sol(p);
                }
            }

            // Set vector and sacalar the solutions
            a = sol.subvec(mat_unknow, size(a));

            scalar = nominal_value;

            double error = arma::norm(X*sol - y, 2) / norm(sol, 2); //  L2 norm

            // return relative error
            return error;


        }

    
        /**
         * @brief estimateError
         * @param cost
         * @param state
         * @param samples
         * @param A
         * @param a
         * @param scalar
         * @return
         */
        double estimateError(std::function<double (const Inputs&)>cost,
                             const Inputs&mean,
                             const mat&samples,
                             InputsMat&A,
                             Inputs&a,
                             double&scalar)
        {
            double error = 0.0;
            Inputs  diff;

            for(unsigned int i=0; i<samples.n_cols; i++)
            {
                Inputs sample   = samples.col(i);
                diff            = sample - mean;

                error += 0.5*(diff.adjoint()*A*diff)[0] + (a.transpose()*diff)[0] + scalar - cost(sample);
            }

            return std::abs(error)/samples.n_cols;
        }

        /**
         * @brief setGaussianSampling
         * @param in
         */
        void setGaussianSampling(const bool in)
        {
            gaussian_sampling = in;
        }

    protected:

        bool gaussian_sampling;
        const unsigned int samples_factor;
        const uint64_t seed;


};


template<int Dim>
void compareRegression(std::function<double (const vec::fixed<Dim>&)>cost,
                       const vec::fixed<Dim>&mean,
                       const mat::fixed<Dim, Dim>&covar,
                       const mat::fixed<Dim, Dim>&M,
                       const vec::fixed<Dim>&m,
                       std::vector<double>&errors,
                       const int&samples_count=100,
                       const double factor = 0.5)
{

    // Check regression
    errors.resize(samples_count);

    const int ncols = covar.n_cols;
    arma::mat Y     = arma::randn(samples_count, ncols);
    mat samples     = arma::repmat(mean, 1, samples_count).t() + Y * arma::chol(covar);

    const double value_mean = cost(mean);

    for(unsigned int i=0; i<samples.n_cols; i++)
    {
        vec::fixed<Dim>vec = samples.col(i);

        const double real   = cost(vec);

        vec                 -= mean;

        errors[i]           = factor*(vec.transpose()*M*vec)[0] +  (vec.transpose()*m)[0] + value_mean;
        errors[i]           -= real;
    }
}

#endif // QUADRATICREGRESSION_HPP
