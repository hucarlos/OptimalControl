#ifndef ILQR_HPP
#define ILQR_HPP

#include <armadillo>

#include <LQR/BasicLQR.hpp>

template<uword xDim, uword uDim>
class iLQR : public BasicLQR<xDim, uDim>
{
    public:
        typedef typename BasicLQR<xDim, uDim>::State State;
        typedef typename BasicLQR<xDim, uDim>::Control Control;

        typedef typename BasicLQR<xDim, uDim>::StateMatrix StateMatrix;
        typedef typename BasicLQR<xDim, uDim>::ControlMatrix ControlMatrix;

        typedef typename BasicLQR<xDim, uDim>::ExtendedState ExtendedState;
        typedef typename BasicLQR<xDim, uDim>::ExtendedStateMatrix ExtendedStateMatrix;

        typedef typename BasicLQR<xDim, uDim>::ControlStateMatrix ControlStateMatrix;
        typedef typename BasicLQR<xDim, uDim>::StateControlMatrix StateControlMatrix;

    public:
        iLQR(const unsigned int ell,
             const BasicSystem<xDim, uDim>*ptr_system,
             const QuadraticCost<xDim> *ptr_init_cost,
             const BasicSystemCost<xDim, uDim>*ptr_system_cost,
             const QuadraticCost<xDim>*ptr_final_cost,
             const bool VIS = true,
             const std::string&name="iLQR"):

            BasicLQR<xDim, uDim>(ell,
                                 ptr_system,
                                 ptr_init_cost,
                                 ptr_system_cost,
                                 ptr_final_cost,
                                 VIS,
                                 name)
        {

        }

        /**
         * @brief estimate
         * @param startState
         * @param maxIter
         * @param delta
         * @param nominalU
         * @return
         */
        double estimate(const State&startState,
                        const unsigned int&maxIter,
                        const double&delta,
                        const std::vector<Control>&nominalU)
        {

            const unsigned int ell = this->getHorizont();

            this->clearAll();

            this->setl(nominalU);

            std::vector< State >   xHat(ell + 1, zeros<vec>(xDim));
            std::vector< Control > uHat(ell, zeros<vec>(uDim));

            double cost     = -std::log(0.0);
            double progress = 0.0;
            double alpha    = 0.0;

            for (unsigned int iter = 0; iter < maxIter; ++iter)
            {
                this->setIteratios(iter);

                // Line search
                cost = this->lineSearch(startState, alpha, cost, progress, xHat, uHat, delta);

                this->setProgress(progress);
                this->setAccum(cost);

                if (this->vis)
                {
                    this->printProgress();
                }

                if (std::abs(progress) < delta)
                {
                    this->nominalState     = xHat;
                    this->nominalControls  = uHat;

                    this->setAccum(cost);

                    return cost;
                }


                // Backward pass
                S         = zeros<mat>(xDim, xDim);
                s         = zeros<vec>(xDim);
                scalar_s  = 0.0;

                this->finalCost->quadratize(xHat.at(ell), S, s, scalar_s);

                // Matrixs for quadratization
                ControlStateMatrix C  = zeros<mat>(uDim, xDim);
                StateMatrix D         = zeros<mat>(xDim, xDim);
                ControlMatrix E       = zeros<mat>(uDim, uDim);
                State d               = zeros<vec>(xDim);
                Control e             = zeros<vec>(uDim);

                for (int t = ell-1; t != -1; --t)
                {

                    // Quadratize cost to go
                    quadratizeCost2Go(xHat, uHat, D, E, C, d, e, t);

                    // Get the policy
                    getPolicy(D, E, C, d, e, t);


                } // Horizont for

            } // Convergence for

            this->nominalState    = xHat;
            this->nominalControls = uHat;

            this->setAccum(cost);

            return cost;
        }

        /**
         * @brief quadratizeCost2Go
         * @param xHat
         * @param uHat
         * @param D
         * @param E
         * @param C
         * @param d
         * @param e
         * @param t
         */
        virtual void quadratizeCost2Go(const std::vector<State>&xHat,
                                       const std::vector<Control>&uHat,
                                       StateMatrix&D,
                                       ControlMatrix&E,
                                       ControlStateMatrix&C,
                                       State&d,
                                       Control&e,
                                       int t)
        {
            const StateMatrix               A = this->system->jacobianState(xHat.at(t), uHat.at(t));
            const StateControlMatrix        B = this->system->jacobianControl(xHat.at(t), uHat.at(t));

            const State c = xHat.at(t+1) - A*xHat.at(t) - B*uHat.at(t);

            ControlStateMatrix P                = zeros<mat>(uDim, xDim);
            StateMatrix Q                       = zeros<mat>(xDim, xDim);
            ControlMatrix R                     = zeros<mat>(uDim, uDim);
            State q                             = zeros<vec>(xDim);
            Control r                           = zeros<vec>(uDim);

            this->systemCost->quadratize(xHat.at(t), uHat.at(t), Q, R, P, q, r);

            //Quadratic terms
            C   = B.t() * S *A + P;
            D   = A.t() * S *A + Q;
            E   = B.t() * S *B + R;

            // Linear terms
            d   = A.t() * (s + S*c) + q;
            e   = B.t() * (s + S*c) + r;

        }


        /**
         * @brief getPolicy
         * @param D
         * @param E
         * @param C
         * @param d
         * @param e
         * @param t
         */
        virtual void getPolicy(const StateMatrix&D,
                               const ControlMatrix&E,
                               const ControlStateMatrix&C,
                               const State&d,
                               const Control&e,
                               const int t)
        {
            this->L.at(t) = -solve(E,C);
            this->l.at(t) = -solve(E,e);

            S = D + C.t() * this->L.at(t);

            s = d + C.t() * this->l.at(t);

            scalar_s = 0.0;
        }


        virtual double lineSearch(const State&initState,
                                  double&alpha,
                                  const double&oldCost,
                                  double&progress,
                                  std::vector<State>&xHat,
                                  std::vector<Control>&uHat,
                                  const double delta=1.0e-4)
        {

            const unsigned int ell = this->getHorizont();

            std::vector< Control > uHatNew(ell, zeros<vec>(uDim));
            std::vector< State > xHatNew(ell + 1, zeros<vec>(xDim));


            double newCost  = 0.0;

            alpha    = 1.0;

            do
            {
                newCost = 0.0;

                // init trajectory
                xHatNew.at(0) = initState;

                for (size_t t = 0; t < ell; ++t)
                {

                    uHatNew.at(t)      = (1.0 - alpha)*uHat.at(t) + this->L.at(t) * (xHatNew.at(t)
                                                                                  - (1.0 - alpha)*xHat.at(t)) + alpha * this->l.at(t);
                    xHatNew.at(t+1)    = this->system->move(xHatNew.at(t), uHatNew.at(t));

                    if(t==0)
                    {
                        newCost += this->initCost->evaluate(xHatNew.at(t));
                    }

                    newCost += this->systemCost->evaluate(xHatNew.at(t), uHatNew.at(t));
                }

                newCost += this->finalCost->evaluate(xHatNew.at(ell));

                alpha *= 0.5;

                progress = (oldCost - newCost) / newCost;

            }
            while (!(newCost < oldCost || std::abs(progress) < delta));

            xHat    = xHatNew;
            uHat    = uHatNew;

            return newCost;

        }

    protected:

        State s;
        StateMatrix S;
        double scalar_s;

};

#endif // ILQR_HPP
