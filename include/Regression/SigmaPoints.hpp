#ifndef SIGMAPOINTS_HPP
#define SIGMAPOINTS_HPP

#include <iostream>
#include <Eigen/Core>
#include <Eigen/Dense>

#include <Utils.hpp>

using namespace std;
using namespace Eigen;

template<int xDim>
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
void sigmaPoints(const Matrix<double,xDim,1>&mean,
                 const Matrix<double, xDim, xDim>&covar,
                 MatrixXd&points,
                 const double&alpha,
                 const double&beta = 2.0,
                 const double&kappa = 0.0
                 )
{
    const unsigned int count_point = (2*xDim) + 1;
    points = Matrix<double, xDim, count_point>::Zero();

    Matrix<double,xDim,xDim> sigma2;

    Eigen::LLT<Eigen::Matrix<double,xDim,xDim> > cholSolver(covar);
    // We can only use the cholesky decomposition if
    // the covariance matrix is symmetric, pos-definite.
    // But a covariance matrix might be pos-semi-definite.
    // In that case, we'll go to an EigenSolver
    if (cholSolver.info()==Eigen::Success)
    {
        // Use cholesky solver
         sigma2 = cholSolver.matrixL();
    }
    else
    {
        throw std::runtime_error("Failed computing the Cholesky decomposition. Use solver instead");
    }

    //Calculate scaling factor for all off-center points
    const double lamda  = (alpha * alpha) * (xDim + kappa) - xDim;
    const double c      = xDim + lamda;

    //    calculate the sigma points; that is,
    //    mu
    //    mu + each column of sigma2 * sqrt(c)
    //    mu - each column of sigma2 * sqrt(c)
    //    Each column of points is one of these.

    points.col(0) = mean;

    for(unsigned int i=1; i<=xDim; i++)
    {
        const Matrix<double,xDim,1> temp = std::sqrt(c) * sigma2.col(i-1);
        points.col(i)       = mean + temp;
        points.col(i+xDim)  = mean - temp;
    }

    //    Calculate weights
//    weights_mean    = Matrix<double, count_point, 1>::Ones();
//    weights_mean[0] = lamda / c;

//    for(unsigned int i=1; i<count_point; i++)
//    {
//        weights_mean[i] = 0.5 / c;
//    }

//    weights_cov     = weights_mean;
//    weights_cov[0]  = lamda / c + (1 - alpha * alpha + beta);
}

//template<int xDim>
///**
// * @brief sigmaPoints
// * @param nominal
// * @param alpha
// * @param points
// */
//void sigmaPoints(const Matrix<double, xDim, 1>&nominal,
//                 const double&alpha,
//                 MatrixXd&points)
//{
//    const  int srDim = 1 << xDim;

//    Matrix<double, srDim, xDim> Ones = Matrix<double, srDim, xDim>::Zero();

//    // Get the necessary samples

//    int lenght = 1 << (xDim-1);

//    for(int c=0; c<xDim; c++)
//    {
//        const int s_count   = 1 << (c+1);

//        for(int s=0, odd=0; s<s_count; s++, odd++)
//        {
//            const int i = s     * lenght;
//            const int f = (s+1) * lenght;

//            for(int r=i; r<f; r++)
//            {
//                if(odd%2 == 0)
//                {

//                    Ones(r,c) = 1.0;
//                }
//                else
//                {
//                    Ones(r,c) = -1.0;
//                }
//            }
//        }

//        lenght = (lenght >> 1);
//    }

//    const int total_samples = (1 << xDim) + 2*xDim + 1;

//    points = MatrixXd::Zero(xDim, total_samples);

//    Matrix<double, xDim, xDim> I = Matrix<double, xDim, xDim>::Identity();


//    //Calculate scaling factor for all off-center points
////    const double kappa  = 0.0;
////    const double lamda  = (alpha * alpha) * (xDim + kappa) - xDim;
//////    const double c      = xDim + lamda;

//    MatrixXd sigma = std::sqrt(alpha) * I;

//    // Get the samples
//    points.col(0) = nominal;
//    for(unsigned int i=1; i<=xDim; i++)
//    {
//        const Matrix<double,xDim,1> temp = sigma.col(i-1);

//        points.col(i)       = nominal + temp;
//        points.col(i+xDim)  = nominal - temp;
//    }


//    // Get the sqrt of the matrix
//    const double factor = 1.0/std::sqrt(xDim);
//    sigma = std::sqrt(alpha) * factor * Ones;

//    const unsigned int init = 2*xDim + 1;

//    for(unsigned int i = init, s=0; s<sigma.rows(); i++, s++)
//    {
//        const Matrix<double,xDim,1> temp    = sigma.row(s);
//        points.col(i)                       = nominal + temp;
//    }

////    printMatrix(Ones);
////    printMatrix(points);

//}

#endif // SIGMAPOINTS_HPP
