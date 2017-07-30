#ifndef FUNCTOR_HPP
#define FUNCTOR_HPP

#include <iostream>
#include <armadillo>

using namespace arma;
using namespace std;

template<int xDim>
class Functor
{
    public:
        typedef arma::vec Input;
        typedef arma::mat ImputsMatrix;

    public:
        Functor()
        {

        }

        virtual double evaluate(const Input&imputs) const = 0;



        virtual void quadratize(const Input&in,
                                ImputsMatrix& Q,
                                Input& q,
                                double&scalar) const
        {
            Q.zeros();
            q.zeros();
            scalar  = 0.0;
        }


        virtual double operator()(const Input&imput) const
        {
            return this->evaluate(imput);
        }
};

#endif // FUNCTOR_HPP


