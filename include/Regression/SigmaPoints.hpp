#ifndef SIGMAPOINTS_HPP
#define SIGMAPOINTS_HPP

#include <iostream>
#include <Utils/Utils.hpp>

using namespace std;
using namespace arma;

//template<int xDim>
/**
 * @brief sigmaPoints
 * @param mean          mean and covariance of a multivariate normal
 * @param covar
 * @param alpha         Spread of the sigma points. Typically 1e-3.
 * @param beta          Used to "incorporate prior knowledge of the distribution of the state".
                        2 is optimal is the state is normally distributed.
 * @param kappa
 * @param points        sigma points
 * @param weights_mean
 * @param weights_cov
 */
//void sigmaPoints(const Matrix<double,xDim,1>&mean,
//                 const Matrix<double, xDim, xDim>&covar,
//                 MatrixXd&points,
//                 const double&alpha,
//                 const double&beta = 2.0,
//                 const double&kappa = 0.0
//                 )
//{
//    const unsigned int count_point = (2*xDim) + 1;
//    points = Matrix<double, xDim, count_point>::Zero();

//    Matrix<double,xDim,xDim> sigma2;

//    Eigen::LLT<Eigen::Matrix<double,xDim,xDim> > cholSolver(covar);
//    // We can only use the cholesky decomposition if
//    // the covariance matrix is symmetric, pos-definite.
//    // But a covariance matrix might be pos-semi-definite.
//    // In that case, we'll go to an EigenSolver
//    if (cholSolver.info()==Eigen::Success)
//    {
//        // Use cholesky solver
//         sigma2 = cholSolver.matrixL();
//    }
//    else
//    {
//        throw std::runtime_error("Failed computing the Cholesky decomposition. Use solver instead");
//    }

//    //Calculate scaling factor for all off-center points
//    const double lamda  = (alpha * alpha) * (xDim + kappa) - xDim;
//    const double c      = xDim + lamda;

//    //    calculate the sigma points; that is,
//    //    mu
//    //    mu + each column of sigma2 * sqrt(c)
//    //    mu - each column of sigma2 * sqrt(c)
//    //    Each column of points is one of these.

//    points.col(0) = mean;

//    for(unsigned int i=1; i<=xDim; i++)
//    {
//        const Matrix<double,xDim,1> temp = std::sqrt(c) * sigma2.col(i-1);
//        points.col(i)       = mean + temp;
//        points.col(i+xDim)  = mean - temp;
//    }

//    //    Calculate weights
////    weights_mean    = Matrix<double, count_point, 1>::Ones();
////    weights_mean[0] = lamda / c;

////    for(unsigned int i=1; i<count_point; i++)
////    {
////        weights_mean[i] = 0.5 / c;
////    }

////    weights_cov     = weights_mean;
////    weights_cov[0]  = lamda / c + (1 - alpha * alpha + beta);
//}

template<int xDim>
/**
 * @brief sigmaPoints
 * @param nominal
 * @param alpha
 * @param points
 */
void sigmaPoints(const vec::fixed<xDim>&mean,
                 const vec::fixed<xDim>&radius,
                 arma::mat&points)
{
    const  int srDim = 1 << xDim;

    mat Ones(srDim, xDim, fill::zeros);

    // Get the necessary diagonal samples
    int lenght = 1 << (xDim-1);

    for(int c=0; c<xDim; c++)
    {
        const int s_count   = 1 << (c+1);

        for(int s=0, odd=0; s<s_count; s++, odd++)
        {
            const int i = s     * lenght;
            const int f = (s+1) * lenght;

            for(int r=i; r<f; r++)
            {
                if(odd%2 == 0)
                {

                    Ones(r,c) = 1.0;
                }
                else
                {
                    Ones(r,c) = -1.0;
                }
            }
        }

        lenght = (lenght >> 1);
    }


    // Configure all points
    const int total_samples = srDim + 2*xDim + 1;

    mat sigma         = eye(xDim, xDim);
    sigma.diag()      = arma::sqrt(radius);

    points            = mat(xDim, total_samples, fill::zeros);


    // Get the samples over every axis
    points.col(0) = mean;
    for(unsigned int i=1; i<=xDim; i++)
    {
        const vec::fixed<xDim> temp = sigma.col(i-1);

        points.col(i)       = mean + temp;
        points.col(i+xDim)  = mean - temp;
    }


    // Get the sqrt of the matrix
    Ones *= 1.0/std::sqrt(xDim);

    const unsigned int init = 2*xDim + 1;

    for(unsigned int i = init, s=0; s< Ones.n_rows; i++, s++)
    {
        const vec::fixed<xDim> temp = arma::sqrt(radius) % Ones.row(s).t();
        points.col(i)               = mean + temp;
    }

}

template<int Dim>
class SigmaPoints
{
    public:

        SigmaPoints(): samplesCount( (1 << Dim) + 2*Dim + 1 )
        {
            const  int srDim = (1 << Dim);

            Ones = zeros<mat>(srDim, Dim);

            // Get the necessary diagonal samples
            int lenght = 1 << (Dim-1);

            for(int c=0; c<Dim; c++)
            {
                const int s_count   = 1 << (c+1);

                for(int s=0, odd=0; s<s_count; s++, odd++)
                {
                    const int i = s     * lenght;
                    const int f = (s+1) * lenght;

                    for(int r=i; r<f; r++)
                    {
                        if(odd%2 == 0)
                        {

                            Ones(r,c) = 1.0;
                        }
                        else
                        {
                            Ones(r,c) = -1.0;
                        }
                    }
                }

                lenght = (lenght >> 1);
            }

        }

        /**
        * @brief estimate
        * @param mean
        * @param radius
        * @param points
        */
        void estimate(const vec::fixed<Dim>&mean,
                      const vec::fixed<Dim>&radius,
                      arma::mat&points)
        {
            // Configure all points
            const int total_samples = this->getSamplesCount();

            mat sigma         = eye(Dim, Dim);
            sigma.diag()      = arma::sqrt(radius);

            points            = mat(Dim, total_samples, fill::zeros);

            // Get the samples over every axis
            points.col(0) = mean;
            for(unsigned int i=1; i<=Dim; i++)
            {
                const vec::fixed<Dim> temp = sigma.col(i-1);

                points.col(i)      = mean + temp;
                points.col(i+Dim)  = mean - temp;
            }


            // Get the sqrt of the matrix
            const double fac = 1.0/std::sqrt(Dim);

            const unsigned int init = 2*Dim + 1;

            for(unsigned int i = init, s=0; s< Ones.n_rows; i++, s++)
            {
                const vec::fixed<Dim> temp  = fac * arma::sqrt(radius) % Ones.row(s).t();
                points.col(i)               = mean + temp;
            }
        }

        /**
        * @brief getSamplesCount
        * @return
        */
        int getSamplesCount() const
        {
            return samplesCount;
        }

    protected:
        const int samplesCount;
        mat Ones;
};

#endif // SIGMAPOINTS_HPP
