#ifndef UTILS_HPP
#define UTILS_HPP



#include <iostream>
#include <fstream>
#include <iomanip> 

#include <vector>
#include <chrono>
#include <utility>

#define ARMA_DONT_PRINT_ERRORS
#include <armadillo>

using namespace arma;

typedef std::chrono::high_resolution_clock::time_point TimeVar;

#define duration(a) std::chrono::duration_cast<std::chrono::milliseconds>(a).count()
#define timeNow() std::chrono::high_resolution_clock::now()

/**
 * @brief regularize
 * @param Q
 * @param epsilon
 */
template <uword aDim>
void regularize0(mat &Q, const double&epsilon)
{


    try
    {
        // Use eigen solver
        vec eigval;
        mat eigvec;

        if( eig_sym(eigval, eigvec, Q, "std"))
        {
            for (unsigned int i = 0; i < aDim; ++i)
            {

                if (eigval(i) < 0.0)
                {
                    eigval(i) = epsilon;
                }

                arma::mat D(aDim, aDim, fill::zeros);
                D.diag() = eigval;

                Q = eigvec * D * eigvec.t();
            }
        }

        else
        {
            throw(std::logic_error("Not eigen decomposition"));
        }

    }
    catch(std::logic_error&e)
    {
        throw std::logic_error(e.what());
    }
}

/**
 * @brief regularize
 * @param Q
 * @param epsilon
 * @param factor
 * @return
 */
template <uword aDim>
void regularize(mat &Q, const double&epsilon, const double&factor)
{


    try
    {
        // Use eigen solver
        vec eigval;
        mat eigvec;


        if( eig_sym(eigval, eigvec, Q, "std"))
        {
            double max = eigval.max();

            if(max <= 0.0)
            {
                for (unsigned int i = 0; i < aDim; ++i)
                {
                   eigval(i) =  eigval(i) - max;

                   if(eigval(i) <= 0.0)
                   {
                       eigval(i) = factor;
                   }
                }

            }

            else
            {
                for (unsigned int i = 0; i < aDim; ++i)
                {

                    if (eigval(i) < epsilon)
                    {
                        eigval(i) = factor;
                    }

                }
            }

            arma::mat D(aDim, aDim, fill::zeros);
            D.diag() = eigval;

            Q = eigvec * D * eigvec.t();
        }

        else
        {
            throw(std::logic_error("Not eigen decomposition regularization"));
        }

    }
    catch(std::logic_error&e)
    {
        throw std::logic_error(e.what());
    }
}


/**
 * @brief regularize2
 * @param Q
 * @param epsilon
 * @param factor
 */
template <uword aDim>
void regularize2(mat &Q, const double&epsilon, const double&factor)
{


    try
    {
        // Use eigen solver
        vec eigval;
        mat eigvec;

        if( eig_sym(eigval, eigvec, Q, "std"))
        {

            for (unsigned int i = 0; i < aDim; ++i)
            {

                if (eigval(i) < epsilon)
                {
                    eigval(i) = factor;
                }

            }

            arma::mat D(aDim, aDim, fill::zeros);
            D.diag() = eigval;

            Q = eigvec * D * eigvec.t();
        }

        else
        {
            throw(std::logic_error("Not eigen decomposition"));
        }

    }
    catch(std::logic_error&e)
    {
        throw std::logic_error(e.what());
    }
}

template <uword aDim>
void regularizeMax(mat &Q, const double&epsilon, const double&factor)
{


    try
    {
        // Use eigen solver
        vec eigval;
        mat eigvec;


        double max = std::numeric_limits<double>::min();

        if( eig_sym(eigval, eigvec, Q, "std"))
        {

            for(unsigned int e=0; e< aDim; e++)
            {
                max = std::max(std::abs(eigval(e)), factor);
            }

            for (unsigned int i = 0; i < aDim; ++i)
            {

                if (eigval(i) < epsilon)
                {
                    eigval(i) = max;
                }

            }

            arma::mat D(aDim, aDim, fill::zeros);
            D.diag() = eigval;

            Q = eigvec * D * eigvec.t();
        }

        else
        {
            throw(std::logic_error("Not eigen decomposition"));
        }

    }
    catch(std::logic_error&e)
    {
        throw std::logic_error(e.what());
    }
}



/**
 * @brief toZero
 * @param M
 * @param epsilon
 */
template<uword xDim, uword yDim>
void toZero(mat&M, const double&epsilon=1.0e-9)
{
    for(unsigned int i=0; i<M.n_rows; i++)
    {
        for(unsigned int j=0; j<M.n_cols; j++)
        {
            if( std::abs(M(i,j)) <= epsilon )
            {
                M(i,j) = 0.0;
            }
        }
    }
}

/**
 * @brief printPath
 * @param path
 * @param filename
 * @return
 */
template<int xDim>
bool printPath(const std::vector< arma::vec::fixed<xDim> >&path,
               const std::string&filename)
{
    FILE * pFile;
    pFile = fopen (filename.c_str(),"w");

    if(pFile == NULL)
    {
        return false;
    }

    for(unsigned int i=0; i<path.size(); i++)
    {
        const arma::vec::fixed<xDim>&x = path[i];

        for(unsigned int j=0; j<xDim; j++)
        {
            fprintf(pFile, "%lf\t", x(j));
        }

        fprintf(pFile, "\n");
    }

    fclose(pFile);

    return true;
}


/**
 * @brief printPathControls
 * @param path
 * @param controls
 * @param filename
 * @return
 */
template<int xDim, int uDim>
bool printPathControls(const std::vector< arma::vec::fixed<xDim> >&path,
                       const std::vector< arma::vec::fixed<uDim> >&controls,
                       const std::string&filename)
{
    FILE * pFile;
    pFile = fopen (filename.c_str(),"w");
    
    if(pFile == NULL)
    {
        return false;
    }
    
    for(unsigned int i=0; i<controls.size(); i++)
    {
        const arma::vec::fixed<xDim>&x = path[i];
        const arma::vec::fixed<uDim>&c = controls[i];

        
        for(int j=0; j<xDim; j++)
        {
            fprintf(pFile, "%lf\t", x(j));
        }
        
        for(int j=0; j<uDim; j++)
        {
            fprintf(pFile, "%lf\t", c(j));
        }
        
        fprintf(pFile, "\n");
    }
    
    
    // Last point
    const int ell = path.size()-1;
    const arma::vec::fixed<xDim>&xell = path[ell];
    
    for(int j=0; j<xDim; j++)
    {
        fprintf(pFile, "%lf\t", xell(j));
    }
    
    for(int j=0; j<uDim; j++)
    {
        fprintf(pFile, "%lf\t", 0.0);
    }
    
    fclose(pFile);
    
    return true;
}


/**
 * @brief printPathControls
 * @param path
 * @param controls
 * @param filename
 * @return
 */
template<int xDim, int uDim>
bool loadPathControls(std::vector< arma::vec::fixed<xDim> >&path,
                      std::vector< arma::vec::fixed<uDim> >&controls,
                       const std::string&filename)
{
    FILE * pFile;
    pFile = fopen (filename.c_str(),"r");

    if(pFile == NULL)
    {
        return false;
    }

    for(unsigned int i=0; i<controls.size(); i++)
    {
        const arma::vec::fixed<xDim>&x = path[i];
        const arma::vec::fixed<uDim>&c = controls[i];


        for(int j=0; j<xDim; j++)
        {
            fscanf(pFile, "%lf\t", x(j));
        }

        for(int j=0; j<uDim; j++)
        {
            fscanf(pFile, "%lf\t", c(j));
        }

        fscanf(pFile, "\n");
    }


    // Last point
    const int ell = path.size()-1;
    const arma::vec::fixed<xDim>&xell = path[ell];

    for(int j=0; j<xDim; j++)
    {
        fscanf(pFile, "%lf\t", xell(j));
    }

    for(int j=0; j<uDim; j++)
    {
        fscanf(pFile, "%lf\t", 0.0);
    }

    fclose(pFile);

    return true;
}

template<uword xDim, uword uDim>
void vector2ArmaMat(const std::vector< arma::vec::fixed<xDim> >&path,
                    const std::vector< arma::vec::fixed<uDim> >&controls,
                    arma::mat&out)
{
    const uword columns = xDim + uDim;
    const uword rows    = path.size();

    out = zeros<mat>(rows, columns);

    for(unsigned int i=0; i<controls.size(); i++)
    {
        const arma::vec::fixed<xDim>&x = path[i];
        const arma::vec::fixed<uDim>&c = controls[i];

        for(int j=0; j<xDim; j++)
        {
            out(i, j) = x(j);
        }

        for(int j=0; j<uDim; j++)
        {
            out(i, j+xDim) = c(j);
        }
    }

    // Last point
    const int ell = path.size()-1;
    const arma::vec::fixed<xDim>&xell = path[ell];

    for(int j=0; j<xDim; j++)
    {
        out(ell, j) = xell(j);
    }

}


/**
 * @brief funcTime
 * @param func
 * @param args
 * @return
 */
template<class F, typename... Args>
double funcTime(F&func,  Args&&... args)
{
    TimeVar t1=timeNow();
    func.estimate(std::forward<Args>(args)...);
    return duration(timeNow()-t1);
}


/**
 From debugger print matrix
 
 @param mat&in input matrix
 */
inline void print_matrix(const mat&in)
{
    for(unsigned int i=0; i<in.n_rows; i++)
    {
        for(unsigned int j=0; j<in.n_cols; j++)
        {
            std::cout << std::left << std::setw(9) << in(i,j) << " ";
        }
        std::cout << std::endl;
    }
    std::cout << std::endl << std::endl;
    
}


/**
 Debbuger print vector
 
 @param vec&in input vector
 */
inline void print_vector(const vec&in)
{
    for(unsigned int i=0; i<in.n_elem; i++)
    {
        std::cout << std::left << std::setw(13) << in(i) << " ";
    }
    std::cout << std::endl << std::endl;
    
    
}

/**
 * @brief operator !
 * @param p
 * @return
 */
template<uword Dim>
mat::fixed<Dim, Dim> operator!(const mat::fixed<Dim, Dim>&p)
{
    const mat::fixed<Dim, Dim>m = arma::solve(p, eye<mat>(Dim, Dim));
    return 0.5 * (m + m.t());
}


/**
 * @brief skewSymmetric
 * @param vector
 * @return
 */
inline mat::fixed<3,3> skewSymmetric(const vec::fixed<3>& vector)
{
    mat::fixed<3,3> result = zeros<mat>(3,3);
    result(0,1) = -vector(2);  result(0,2) =  vector(1);
    result(1,0) =  vector(2);  result(1,2) = -vector(0);
    result(2,0) = -vector(1);  result(2,1) =  vector(0);

    return result;
}

template <uword Dim>
inline mat::fixed<Dim, Dim> expSkewSym(const vec::fixed<Dim>&q)
{
    const mat::fixed<Dim, Dim>skew = skewSymmetric(q);
    const double x = arma::norm(q, 2);

    double sine     = (std::sin(x)/x);
    double cosine   = (1.0 - std::cos(x))/(x*x);

    if(x <= 1.0e-6)
    {
            sine    = 1.0;
            cosine  = 0.5;
    }

    mat::fixed<Dim, Dim> result = eye<mat>(Dim, Dim) + sine*skew + cosine*skew*skew;
    return result;
}

// Matrix exponentiation
#define _MATRIX_B0 1729728e1
#define _MATRIX_B1 864864e1
#define _MATRIX_B2 199584e1
#define _MATRIX_B3 2772e2
#define _MATRIX_B4 252e2
#define _MATRIX_B5 1512e0
#define _MATRIX_B6 56e0
#define _MATRIX_B7 1e0
#define _NORMLIM 9.504178996162932e-1

template <uword Dim>
inline mat::fixed<Dim, Dim> exp(const mat::fixed<Dim, Dim>& q)
{
    mat::fixed<Dim, Dim> A(q);
    int s = (int) std::max(double(0), ceil(log(norm(A)/_NORMLIM)*M_LOG2E));

    const mat::fixed<Dim, Dim>identity = eye<mat>(Dim,Dim);

    A /= pow(2.0,s);
    mat::fixed<Dim, Dim> A2(A*A);
    mat::fixed<Dim, Dim> A4(A2*A2);
    mat::fixed<Dim, Dim> A6(A2*A4);
    mat::fixed<Dim, Dim> U( A*(A6*_MATRIX_B7 + A4*_MATRIX_B5 + A2*_MATRIX_B3 + identity*_MATRIX_B1) );
    mat::fixed<Dim, Dim> V( A6*_MATRIX_B6 + A4*_MATRIX_B4 + A2*_MATRIX_B2 + identity*_MATRIX_B0 );
    mat::fixed<Dim, Dim> R7 = arma::solve((V - U),(V + U));

    for (int i = 0; i < s; ++i)
    {
        R7 = R7*R7;
    }
    return R7;
}


template<uword Dim>
inline void vec2DiagMat(const vec::fixed<Dim>&in, mat::fixed<Dim, Dim>&out)
{
    for(unsigned int i=0; i<Dim; i++)
    {
        for(unsigned int j=0; j<Dim; j++)
        {
            out(i,j) = (i==j)?(in(i)):(0.0);
        }
    }
}



#endif // UTILS_HPP
