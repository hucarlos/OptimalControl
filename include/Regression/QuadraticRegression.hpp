#ifndef QUADRATICREGRESSION_HPP
#define QUADRATICREGRESSION_HPP

#include <iostream>
#include <armadillo>

#include <Regression/SigmaPoints.hpp>

#include "tbb/parallel_for.h"
#include "tbb/task_scheduler_init.h"

using namespace std;
using namespace arma;

template<uword Dim>
class QuadraticRegression
{
    public:
        typedef arma::vec::fixed<Dim>Inputs;
        typedef arma::mat::fixed<Dim, Dim>InputsMat;

        enum SAMPLING_MODE {GAUSSIAN_S, SIGMA_S, ELLIPSOID_S};

        /**
         * @brief QuadraticRegression
         */
        QuadraticRegression(): samplingMode(GAUSSIAN_S), samplingFactor(1),
            seed(0),
            _unknown(((Dim * (Dim+1))/2) + Dim + 1)
        {
            this->updateSamplesMat();
        }

        /**
         * @brief QuadraticRegression
         * @param mode
         */
        QuadraticRegression(SAMPLING_MODE mode): samplingMode(mode),
            samplingFactor(1),
            seed(0),
            _unknown(((Dim * (Dim+1))/2) + Dim + 1)
        {
            this->updateSamplesMat();
        }

        /**
         * @brief QuadraticRegression
         * @param mode
         * @param sf
         */
        QuadraticRegression(SAMPLING_MODE mode, const unsigned int&sf):  samplingMode(mode), samplingFactor(sf),
            seed(0),
            _unknown(((Dim * (Dim+1))/2) + Dim + 1)
        {
          this->updateSamplesMat();
        }

        /**
         * @brief QuadraticRegression
         * @param mode
         * @param s
         * @param sf
         */
        QuadraticRegression(SAMPLING_MODE mode,  const unsigned int&sf, const unsigned int&s):
            samplingMode(mode), samplingFactor(sf),
            seed(s),
            _unknown(((Dim * (Dim+1))/2) + Dim + 1)
        {
            this->updateSamplesMat();
        }

    public:

        /**
         * @brief getRegression
         * @param cost
         * @param mean
         * @param A
         * @param a
         * @param scalar
         * @param radius
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
            arma::mat samples;
            InputsMat cov;
            vec2DiagMat(radius, cov);


            if(samplingMode == GAUSSIAN_S)
            {

                samples = mat(Dim, samplesCount, fill::zeros);
                samples = (arma::repmat(mean, 1, samplesCount).t() + Y * arma::chol(cov)).t();
            }

            else if(samplingMode == ELLIPSOID_S)
            {

                mat esamples;
                if(_parallel)
                {
                    ellipPoints.getSamples(mean, cov, esamples);
                }
                else
                {
                    ellipPoints.getSamples(mean, cov, esamples);
                }
                samples = esamples.t();
            }

            else if(samplingMode == SIGMA_S)
            {

                Inputs startRadius = radius;
                const int count_points   = std::max(int(std::ceil((Dim + 3)/8)), 1);
                const int samples4Sigma  =  (1 << Dim) + 2*Dim;

                const unsigned int count_samples  = samplingFactor * count_points;

                arma::mat points;

                samples = mat(Dim, samples4Sigma*count_samples, fill::zeros);


                for(unsigned int i=0; i<count_samples; i++)
                {
                    //sigmaPoints<xDim>(mean, radius, points);

                    sPoints.estimate(mean, radius, points);

                    samples.submat(0, i*samples4Sigma, size(Dim, samples4Sigma)) = points.submat(0,1,size(Dim, samples4Sigma));

                    startRadius *= 0.5;
                }

            }

            double relative_error = this->solve(cost, mean, samples, A, a, scalar, lambda);
            return relative_error;
        }

        /**
         * @brief getRegression
         * @param cost
         * @param mean
         * @param A
         * @param a
         * @param scalar
         * @param alpha
         * @param lambda
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
            Inputs radius = alpha * ones<vec>(Dim);
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


            const int mat_unknow    = (Dim*(Dim + 1)/2);
            const int unknow        =  mat_unknow + Dim;

            if(unknow <= 0)
            {
                throw std::runtime_error("Bad dimensions for unknow on solve");
            }

            const int samples_count = samples.n_cols;


            mat X(samples_count, unknow, fill::zeros);
            vec y(samples_count, fill::zeros);


            // Cicle for every sample
            const double nominal_value = cost(mean);

            if(_parallel)
                buildParallel(cost, mean, nominal_value, samples, X, y);
            else
                buildSerial(cost, mean, nominal_value, samples, X, y);


            // Obtain the solution
            const mat::fixed<unknow, unknow>I = eye<mat>(unknow, unknow);

            mat::fixed<unknow, unknow> L    = X.t() * X;
            vec::fixed<unknow>b             = X.t() * y;

            mat::fixed<unknow, unknow> M    = (L + lambda*I);
            vec::fixed<unknow> sol;

            mat R;
            if(arma::chol(R, M))
            {
                // Solve the system
                vec y   = arma::solve(trimatl(R.t()), b);
                sol     = arma::solve(trimatu(R), y);
            }
            else
            {
                sol = arma::solve(M, b, solve_opts::equilibrate);
            }


            // Put results on A
            for(unsigned int i=0, p=0; i<Dim; i++)
            {
                for(unsigned int j=i; j<Dim; j++, p++)
                {
                    A(i,j) = A(j,i) = sol(p);
                }
            }

            // Set vector and sacalar the solutions
            a = sol.subvec(mat_unknow, arma::size(a));

            scalar = nominal_value;

            double error = arma::norm(X*sol - y, 2) / norm(y, 2); //  L2 norm

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
        void setSamplingMode(SAMPLING_MODE mode)
        {
            samplingMode = mode;
        }

        /**
         * @brief setSamplingFactor
         * @param sf
         */
        void setSamplingFactor(const unsigned int&sf)
        {
            samplingFactor = sf;
            this->updateSamplesMat();

        }

        /**
         * @brief setSeed
         * @param s
         */
        void setSeed(const unsigned int&s)
        {
            seed = s;
        }

        /**
         * @brief setParallel
         * @param in
         */
        void setParallel(const bool&in)
        {
            _parallel = in;
        }

        /**
         * @brief parallel
         * @return
         */
        bool parallel() const
        {
            return _parallel;
        }


    protected:

        /**
         * @brief buildSerial
         * @param cost
         * @param mean
         * @param nominal_value
         * @param samples
         * @param X
         * @param y
         */
        void buildSerial(std::function<double (const Inputs&)>cost, const vec&mean,
                         const double&nominal_value, const mat&samples, mat&X, vec&y)
        {
            for(unsigned int i=0; i<samples.n_cols; i++)
            {
                const Inputs sample = samples.col(i);
                const double value  = cost(sample);
                const Inputs diff   = sample - mean;

                // Cicle for all the matrix values
                int p = 0;
                for(unsigned int r=0; r<Dim; r++)
                {
                    for(unsigned int c=r; c<Dim; c++, p++)
                    {
                        X(i,p)  = diff(r)*diff(c);

                        X(i,p) *= 0.5;

                        if(c>r)
                        {
                            X(i,p) *= 2.0;
                        }
                    }
                }

                for(unsigned int r=0; r<Dim; r++, p++)
                {
                    X(i,p) = diff(r);
                }

                y(i) = (value - nominal_value);

            }
        }


        /**
         * @brief buildParallel
         * @param cost
         * @param mean
         * @param nominal_value
         * @param samples
         * @param X
         * @param y
         */
        void buildParallel(std::function<double (const Inputs&)>cost, const vec&mean,
                           const double&nominal_value, const mat&samples, mat&X, vec&y)
        {
            tbb::task_scheduler_init init(tbb::task_scheduler_init::automatic);
            tbb::parallel_for(
                        tbb::blocked_range<size_t>(0, samples.n_cols),
                        [&X,&y,&mean,&nominal_value, &samples, &cost](const tbb::blocked_range<size_t>& r)
            {
                for (size_t i=r.begin();i<r.end();++i)
                {
                    const Inputs sample = samples.col(i);
                    const double value  = cost(sample);
                    const Inputs diff   = sample - mean;

                    // Cicle for all the matrix values
                    int p = 0;
                    for(unsigned int r=0; r<Dim; r++)
                    {
                        for(unsigned int c=r; c<Dim; c++, p++)
                        {
                            X(i,p)  = diff(r)*diff(c);

                            X(i,p) *= 0.5;

                            if(c>r)
                            {
                                X(i,p) *= 2.0;
                            }
                        }
                    }

                    for(unsigned int r=0; r<Dim; r++, p++)
                    {
                        X(i,p) = diff(r);
                    }

                    y(i) = (value - nominal_value);

                }

            });
        }

        /**
         * @brief updateSamplesMat
         */
        void updateSamplesMat()
        {
            samplesCount = samplingFactor * _unknown;

            // For Gaussian mode
            arma_rng::set_seed(seed);
            Y   = arma::randn(samplesCount, Dim);

            ellipPoints.setSamplesCount(samplesCount);

        }



    protected:

        SAMPLING_MODE samplingMode;
        unsigned int samplingFactor;
        unsigned int seed;

        bool _parallel;



        const unsigned int _unknown;
        unsigned int samplesCount;

        arma::mat Y;

        EllipsoidPoints<Dim> ellipPoints;

        SigmaPoints<Dim>sPoints;


};


/**
 * @brief compareRegression
 * @param cost
 * @param mean
 * @param covar
 * @param M
 * @param m
 * @param errors
 * @param samples_count
 * @param factor
 */
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
