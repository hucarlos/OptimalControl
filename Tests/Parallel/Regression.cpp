#include "tbb/parallel_for.h"
#include "tbb/task_scheduler_init.h"
#include <iostream>
#include <vector>

#include <armadillo>
using namespace arma;

#define DIM 5

arma::vec cost(const vec&input)
{
    return 2*input;
}

int main(int,char**)
{
    
    //tbb::task_scheduler_init init;  // Automatic number of threads
    tbb::task_scheduler_init init(tbb::task_scheduler_init::default_num_threads());  // Explicit number of threads

    const int samples_count = 300;

    vec mean(DIM, fill::randu);
    mat samples(DIM, samples_count, fill::randn);

    //Serial function

    std::vector<double>norm1(samples_count, 0.0);
    std::vector<double>norm2(samples_count, 0.0);


    for(int i=0; i<samples_count; i++)
    {
        const vec sample = samples.col(i);
        const vec value  = cost(sample);
        norm1.at(i)      = norm(value - mean, 2);
    }

    //Parallel version

    tbb::parallel_for(
                      tbb::blocked_range<size_t>(0,samples_count),
                      [&norm2,&samples, &mean](const tbb::blocked_range<size_t>& r) {
                          for (size_t i=r.begin();i<r.end();++i)
                          {
                              const vec sample = samples.col(i);
                              const vec value  = cost(sample);
                              norm2.at(i)      = norm(value - mean, 2);
                          }
                      }
                      );
    
    for(int i=0; i<samples_count; i++)
    {
        std::cerr << norm1[i]-norm2[i]<< std::endl;
    }
    
    return 0;
}
