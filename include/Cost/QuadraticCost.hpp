#ifndef QUADRATICCOST_HPP
#define QUADRATICCOST_HPP

#include <iostream>

#include <Cost/Functor.hpp>

using namespace arma;

template<std::size_t xDim>
class QuadraticCost : public Functor<xDim>
{
    public:
        typedef typename Functor<xDim>::Input Input;

        typedef typename Functor<xDim>::ImputsMatrix ImputsMatrix;

    public:
        QuadraticCost(): _nominal(zeros<vec>(xDim)), _Q(eye<mat>(xDim, xDim))
        {


        }

        QuadraticCost(const Input&nominal): _nominal(nominal), _Q(eye<mat>(xDim, xDim))
        {

        }

        QuadraticCost(const ImputsMatrix&Q): _nominal(Input::Zero()), _Q(Q)
        {

        }

        QuadraticCost(const Input&nominal, const ImputsMatrix&Q): _nominal(nominal), _Q(Q)
        {

        }

        QuadraticCost(const double&lambda): _nominal(zeros<vec>(xDim)),  _Q(lambda * eye<mat>(xDim, xDim))
        {

        }

        QuadraticCost(const QuadraticCost<xDim>&other): _nominal(other._nominal), _Q(other._Q)
        {

        }

        double evaluate(const Input&input) const
        {
            const Input diff = (input - _nominal);
            return 0.5*(dot(diff,  (_Q * diff)));
        }

        Input gradient(const Input&input) const
        {
            const Input diff = (input - _nominal);
            return _Q * diff;
        }

        ImputsMatrix hessian() const
        {
            return _Q;
        }



        /**
         * @brief quadratize
         * @param inputs
         * @param Q
         * @param q
         */
        void quadratize(const Input&inputs, ImputsMatrix&Q, Input&q, double&scalar) const
        {
            Q       += hessian();
            q       += -(_Q * _nominal);
            scalar  += 0.5*(dot(_nominal,  (_Q * _nominal)));
        }

        /**
         * @brief constant
         * @return
         */
        double constant() const
        {
            return 0.5*(dot(_nominal, (_Q * _nominal)));
        }

    private:

        const Input _nominal;
        const ImputsMatrix _Q;

};

#endif // QUADRATICCOST_HPP
