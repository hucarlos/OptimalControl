#ifndef BASICLQR_HPP
#define BASICLQR_HPP

#include <iostream>
#include <iomanip>
#include <armadillo>

#include <Cost/BasicSystemCost.hpp>
#include <Cost/QuadraticCost.hpp>
#include <System/BasicSystem.hpp>

template<uword xDim, uword uDim>
class BasicLQR
{
    public:
        typedef arma::vec::fixed< xDim >State;
        typedef arma::vec::fixed< uDim >Control;

        typedef arma::mat::fixed< xDim, xDim >StateMatrix;
        typedef arma::mat::fixed< uDim, uDim >ControlMatrix;

        typedef arma::vec::fixed< xDim + uDim>ExtendedState;
        typedef arma::mat::fixed< xDim + uDim, xDim + uDim>ExtendedStateMatrix;

        typedef arma::mat::fixed< uDim, xDim>ControlStateMatrix;
        typedef arma::mat::fixed< xDim, uDim>StateControlMatrix;

    public:

        BasicLQR(const unsigned int ell,
                   const BasicSystem<xDim, uDim>*ptr_system,
                   const QuadraticCost<xDim> *ptr_init_cost,
                   const BasicSystemCost<xDim, uDim>*ptr_system_cost,
                   const QuadraticCost<xDim>*ptr_final_cost,
                   const bool VIS = true,
                   const std::string&name="BasicLQR"):

            horizont(ell),
            system(ptr_system),
            initCost(ptr_init_cost),
            systemCost(ptr_system_cost),
            finalCost(ptr_final_cost),
            vis(VIS),
            className(name),

            L(std::vector<ControlStateMatrix>(ell)),
            l(std::vector<Control>(ell)),

            nominalState(std::vector<State>(ell+1)),
            nominalControls(std::vector<Control>(ell))

        {

        }

        ~BasicLQR()
        {
            system      = NULL;
            initCost    = NULL;
            finalCost   = NULL;
            systemCost  = NULL;
        }


        virtual double estimate(const State&startState,
                                const unsigned int&maxIter,
                                const double&delta,
                                const std::vector<Control>&nominalU) = 0;



        /**
         * @brief estimatePath
         * @param xStart
         * @param path
         * @param controls
         * @return
         */
        double estimatePath(const State&xStart)
        {
            const unsigned int ell = getHorizont();

            std::fill(nominalState.begin(),     nominalState.end(),     zeros<vec>(xDim));
            std::fill(nominalControls.begin(),  nominalControls.end(),  zeros<vec>(uDim));

            State x             = xStart;
            nominalState.at(0)  = x;

            double cost = 0.0;

            cost += initCost->evaluate(x);

            for (unsigned int t = 0; t < ell; ++t)
            {
                Control u       = L.at(t)*x + l.at(t);
                nominalControls.at(t)  = u;

                cost += systemCost->evaluate(x, u);

                x = system->move(x, u);
                nominalState.at(t+1) = x;

            }

            cost += finalCost->evaluate(x);
            nominalState.at(ell) = x;

            return cost;
        }

        /**
         * @brief clearAll
         */
        void clearAll()
        {
            std::fill(L.begin(), L.end(), zeros<mat>(uDim, xDim));
            std::fill(l.begin(), l.end(), zeros<vec>(uDim));

            std::fill(nominalState.begin(), nominalState.end(), zeros<vec>(xDim));

            std::fill(nominalControls.begin(), nominalControls.end(), zeros<vec>(uDim));

            setAccum(0.0);
            setIteratios(0);
            setTime(0);
            setProgress(0.0);
        }

        void setl(const std::vector<Control>&in)
        {
            if(in.size() == getHorizont())
            {
                l = in;
            }
            else
            {
                throw std::runtime_error("Bad open loop control size");
            }
        }

    public:

        void printProgress()
        {
            std::cout <<"Iter: "        << std::left << std::setw(8)  << std::setprecision(7) << iterations() <<"  ";
            std::cout <<"Progress: "    << std::left << std::setw(15) << std::setprecision(7) << getProgress() <<"  ";
            std::cout <<"Cost: "        << std::left << std::setw(15) << std::setprecision(7) << getAccum()  <<endl;
        }
    
        virtual void  printConfiguration()
        {
            
        }

        /**
         * @brief getIters
         * @return
         */
        unsigned int iterations() const
        {
            return iters;
        }

        /**
         * @brief getFinalCost
         * @return
         */
        double getAccum() const
        {
            return accum;
        }

        /**
         * @brief BasicLQR::getHorizont
         * @return
         */
        unsigned int getHorizont() const
        {
            return horizont;
        }

        double getProgress() const
        {
            return progress;
        }


        void getNominalState(std::vector<State>&out) const
        {
            out = nominalState;
        }

        void getNominalControl(std::vector<Control>&out) const
        {
            out = nominalControls;
        }

        std::string getName() const
        {
            return className;
        }

    protected:

        const unsigned int horizont;

        const BasicSystem<xDim, uDim>*system;
        const QuadraticCost<xDim>*initCost;
        const BasicSystemCost<xDim, uDim>*systemCost;
        const QuadraticCost<xDim>*finalCost;

        const bool vis;
        const std::string className;


        std::vector< ControlStateMatrix  >L;
        std::vector< Control >l;

        std::vector<State>nominalState;
        std::vector<Control>nominalControls;


    protected:

        void setTime(const unsigned int&t)
        {
            time = t;
        }

        void setIteratios(const unsigned int&value)
        {
            iters = value;
        }

        void setAccum(const double&value)
        {
            accum = value;
        }

        void setProgress(const double&value)
        {
            progress = value;
        }



        unsigned int time;
        unsigned int iters;
        double accum;
        double progress;

};

#endif // BASICLQR_HPP






