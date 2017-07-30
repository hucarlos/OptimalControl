#ifndef UTILS_HPP
#define UTILS_HPP



#include <iostream>
#include <fstream>
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
 * @param factor
 * @return
 */
template <uword aDim>
int regularize(mat &Q, const double&epsilon, const double&factor)
{

    int counter = 0;

    // Use eigen solver
    vec eigval;
    mat eigvec;

    eig_sym(eigval, eigvec, Q, "std");

    for (int i = 0; i < aDim; ++i)
    {
        if (eigval(i) < epsilon)
        {
            eigval(i) = factor;
            counter++;
        }
    }

    arma::mat D(aDim, aDim, fill::zeros);
    D.diag() = eigval;

    Q = eigvec * D * eigvec.t();

    return counter;
}


/**
 * @brief toZero
 * @param M
 * @param epsilon
 */
template<uword xDim, uword yDim>
void toZero(arma::mat::fixed<xDim, yDim>&M, const double&epsilon=1.0e-06)
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

#endif // UTILS_HPP
