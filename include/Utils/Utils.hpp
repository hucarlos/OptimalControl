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
void regularize(mat &Q, const double&epsilon, const double&factor)
{


    try
    {
        // Use eigen solver
        vec eigval;
        mat eigvec;

        if( eig_sym(eigval, eigvec, Q, "std"))
        {

            for (int i = 0; i < aDim; ++i)
            {
                if(eigval(i)<0)
                {
                    eigval(i) *= -1.0;
                }

                if (eigval(i) < epsilon)
                {
                    eigval(i) = factor;
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


template <uword aDim>
/**
 * @brief regularize
 * @param Q
 * @param epsilon
 * @param factor
 */
void regularize(const mat &Q, mat&Qp, const double&epsilon, const double&factor)
{


    try
    {
        // Use eigen solver
        vec eigval;
        mat eigvec;

        if( eig_sym(eigval, eigvec, Q, "std"))
        {

            for (int i = 0; i < aDim; ++i)
            {
//                if(eigval(i)<0)
//                {
//                    eigval(i) *= -1.0;
//                }

                if (eigval(i) < epsilon)
                {
                    eigval(i) = factor;
                }

                arma::mat D(aDim, aDim, fill::zeros);
                D.diag() = eigval;

                Qp = eigvec * D * eigvec.t();
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


#include <array>
#include <vector>
#include <algorithm>

template <typename T, size_t N>
/**
 * @brief generate_rectified_hypercube
 * @return
 */
std::vector<std::array<T, N>> generate_rectified_hypercube()
{
    std::vector<std::array<T, N>> vertices;  // result data
    std::array<T, N> set; // set to permute upon
    // Initialize set to { 0, 1, ..., 1 }
    set[0] = T(0);
    for (size_t i = 1; i < N; ++i) set[i] = T(1);

    // Generate vertices

    // Add all possible permutations of initial set (unrolled loop for i=0)
    do {
        vertices.push_back(set);
    } while (std::next_permutation(set.begin(), set.end()));

    for (size_t i = 1; i < N; ++i) { // do the rest of them

        // Modify set to be seed for next group of permutations:
        // "abusing" the known state after last std::next_permutation call
        set[i-1] = T(-1);    // Here was zero previously
        set[i] = T(0);       // Here was first 1

        // Add all possible permutations of current set
        do {
            vertices.push_back(set);
        } while (std::next_permutation(set.begin(), set.end()));

    }
    return vertices;
}

#endif // UTILS_HPP
