#ifndef BASICSYSTEMCOST_HPP
#define BASICSYSTEMCOST_HPP

#include <iostream>
#include <armadillo>

using namespace arma;


template<uword xDim, uword uDim>
class BasicSystemCost
{
    public:

        typedef arma::vec::fixed<xDim>State;
        typedef arma::vec::fixed<uDim>Control;

        typedef arma::mat::fixed< xDim, xDim>StateMatrix;
        typedef arma::mat::fixed< uDim, uDim>ControlMatrix;

        typedef arma::mat::fixed< uDim, xDim>ControlStateMatrix;

    public:
        BasicSystemCost()
        {

        }

        virtual double evaluate(const State&state, const Control&control) const = 0;

        /**
         * @brief quadratize
         * @param state
         * @param control
         * @param Qt
         * @param Rt
         * @param Pt
         * @param qt
         * @param rt
         * @return
         */
        virtual void quadratize(const State&state, const Control&control,
                                  StateMatrix&Qt, ControlMatrix&Rt, ControlStateMatrix&Pt,
                                  State&qt, Control&rt) const = 0;


        /**
         * @brief operator ()
         * @param input
         * @return
         */
        virtual double operator ()(const vec::fixed<(xDim+uDim)>&input)
        {
            const int indexS = xDim -1;
            const int indexC = (xDim + uDim) -1;
            return this->evaluate(input.subvec(0, indexS), input.subvec(xDim, indexC));
        }

        /**
         * @brief operator ()
         * @param state
         * @param control
         * @return
         */
        virtual double operator()(const State&state, const Control&control) const
        {
            return this->evaluate(state, control);
        }

        /**
         * @brief operator ()
         * @param control
         * @return
         */
        virtual double evaluateControl(const Control&control) const
        {
            return this->evaluate(zeros<vec>(xDim), control);
        }

        /**
         * @brief operator ()
         * @param state
         * @return
         */
        virtual double evaluateState(const State&state) const
        {
            return this->evaluate(state, zeros<vec>(uDim));
        }


    protected:


};

#endif // BASICSYSTEMCOST_HPP

