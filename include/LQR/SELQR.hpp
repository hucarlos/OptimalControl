#ifndef SELQR_HPP
#define SELQR_HPP

#include <iostream>
#include <LQR/BasicLQR.hpp>

template<uword xDim, uword uDim>
class SELQR : public BasicLQR<xDim, uDim>
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

        SELQR(const unsigned int ell,
              const BasicSystem<xDim, uDim>*ptr_system,
              const QuadraticCost<xDim> *ptr_init_cost,
              const BasicSystemCost<xDim, uDim>*ptr_system_cost,
              const QuadraticCost<xDim>*ptr_final_cost,
              const bool VIS = true,
              const std::string&name="SELQR"):

            BasicLQR<xDim, uDim>(ell,
                                 ptr_system,
                                 ptr_init_cost,
                                 ptr_system_cost,
                                 ptr_final_cost,
                                 VIS,
                                 name),

            SBar(std::vector<StateMatrix>(ell+1)),
            sBar(std::vector<State>(ell+1)),
            scalar_sBar(std::vector<double>(ell+1)),

            S(std::vector<StateMatrix>(ell+1)),
            s(std::vector<State>(ell+1)),
            scalar_s(std::vector<double>(ell+1))


        {

        }

    public:
        double estimate(const State&startState,
                        const unsigned int&maxIter,
                        const double&delta,
                        const std::vector<Control>&nominalU)
        {
            // Initialization
            const unsigned int ell = this->getHorizont();


            State xHat          = startState;

            reset(nominalU);

            double newCost  = 0.0;
            double oldCost  = -std::log(0.0);


            for (unsigned int iter = 0; iter < maxIter; ++iter)
            {

                try
                {
                    this->setIteratios(iter);

                    SBar.at(0).zeros();
                    sBar.at(0).zeros();
                    scalar_sBar.at(0)   = 0.0;

                    // ================================ FORWARD PASS ========================================
                    for (unsigned int t = 0; t < ell; ++ t)
                    {

                        this->setTime(t);

                        // Get the nominal point
                        const Control uHat       = this->L.at(t)*xHat + this->l.at(t);

                        // ====================== Quadratize cost-to-come ==================================

                        StateMatrix DBar    = zeros<mat>(xDim, xDim);
                        ControlMatrix EBar  = zeros<mat>(uDim, uDim);

                        ControlStateMatrix CBar = zeros<mat>(uDim, xDim);

                        State dBar      = zeros<vec>(xDim);
                        Control eBar    = zeros<vec>(uDim);

                        double scalar_f = 0.0;

                        quadratizeCost2Come(xHat, uHat, t, DBar, EBar, CBar, dBar, eBar, scalar_f);

                        getForwardPolicy(DBar, EBar, CBar, dBar, eBar, scalar_f, t);

                        // Update nominal state t+1
                        const int update = t+1;
                        estimateNominalState(update, xHat);

                    }

                    // ================================ BACKWARD PASS ========================================

                    S.at(ell).zeros();
                    s.at(ell).zeros();
                    scalar_s.at(ell)  = 0.0;

                    this->finalCost->quadratize(xHat, S.at(ell), s.at(ell), scalar_s.at(ell));

                    // Update nominal state in ell
                    estimateNominalState(ell, xHat);

                    for (int t = ell - 1; t != -1; --t)
                    {
                        this->setTime(t);

                        // Get the nominal point
                        const Control uHat       = this->L.at(t)*xHat + this->l.at(t);

                        // ====================== Quadratize cost-to-come =====================

                        StateMatrix D   = zeros<mat>(xDim,xDim);
                        ControlMatrix E = zeros<mat>(uDim, uDim);

                        ControlStateMatrix C = zeros<mat>(uDim, xDim);

                        State d     = zeros<vec>(xDim);
                        Control e   = zeros<vec>(uDim);

                        double scalar_b = 0.0;

                        quadratizeCost2Go(xHat, uHat, t, D, E, C, d, e, scalar_b);

                        // ====================== Get policy ==================================

                        getBackwardPolicy(D, E, C, d, e, scalar_b, t);

                        // Update nominal state in t
                        estimateNominalState(t, xHat);

                    }

                    // compute cost
                    newCost = this->estimatePath(xHat);

                    this->progress = ((oldCost - newCost) / newCost);

                    this->setAccum(newCost);

                    if(this->vis)
                    {
                        this->printProgress();
                    }


                    // Stop condition
                    if (std::abs(this->getProgress()) <= delta)
                    {
                        ++iter;
                        this->setIteratios(iter);
                        return newCost;
                    }

                    oldCost = newCost;

                    // Set the parameters to the new values for every iteration
                    setParameters();

                    // If something goes wrong try to start again
                    if(std::isnan(this->getAccum()) || (!isfinite(this->getAccum())) /*|| (this->getProgress() == -1.0)*/)
                    {
                        newCost  = 0.0;
                        oldCost  = -std::log(0.0);

                        xHat       = startState;

                        setInitialConditions(nominalU);
                    }
                }
                catch(std::logic_error&e)
                {
                    newCost  = 0.0;
                    oldCost  = -std::log(0.0);

                    xHat       = startState;

                    setInitialConditions(nominalU);
                    continue;
                }

            }

            this->setIteratios(maxIter);
            newCost = this->estimatePath(xHat);
            if(std::isnan(newCost))
            {
                this->setAccum(-std::log(0.0));
            }
            else
            {
                this->setAccum(newCost);
            }
            return newCost;

        }

        /**
         * @brief quadratizeCost2Come
         * @param xHat
         * @param uHat
         * @param t
         * @param DBar
         * @param EBar
         * @param CBar
         * @param dBar
         * @param eBar
         * @param scalar_f
         */
        virtual void quadratizeCost2Come(const State&xHat, const Control&uHat, const int&t,
                                         StateMatrix&DBar,
                                         ControlMatrix&EBar,
                                         ControlStateMatrix&CBar,
                                         State&dBar,
                                         Control&eBar,
                                         double&scalar_f)
        {
            // ====================== Linearize forward movement model ==================================

            StateMatrix ABar          = zeros<mat>(xDim, xDim);
            StateControlMatrix BBar   = zeros<mat>(xDim, uDim);
            State cBar                = zeros<vec>(xDim);

            linearizeForward(xHat, uHat, ABar, BBar, cBar);

            // ====================== Quadratize instant cost  ==================================

            ControlStateMatrix P    = zeros<mat>(uDim, xDim);
            StateMatrix Q           = zeros<mat>(xDim, xDim);
            ControlMatrix R         = zeros<mat>(uDim, uDim);

            State q                 = zeros<vec>(xDim);
            Control r               = zeros<vec>(uDim);

            quadratizeInstantCost(xHat, uHat, t, Q, R, P, q, r);


            // Quadratization

            const StateMatrix SBarQ     = SBar.at(t) + Q;  // SBarQ Symetric
            const State sBarqSBarQcBar  = sBar.at(t) + q + (SBarQ*cBar);

            CBar    = (BBar.t() * SBarQ * ABar) + P*ABar;
            DBar    = (ABar.t() * SBarQ * ABar); // Symetric
            EBar    = (BBar.t() * SBarQ * BBar) + R + (P*BBar) +  (P*BBar).t() ; // Symetric

            // Linear terms
            dBar    = ABar.t() * sBarqSBarQcBar;
            eBar    = BBar.t() * sBarqSBarQcBar + r + P*cBar;

            scalar_f = 0.0; //Not necessary in SELQR
        }

        /**
         * @brief quadratizeCost2Go
         * @param xHat
         * @param uHat
         * @param t
         * @param D
         * @param E
         * @param C
         * @param d
         * @param e
         * @param scalar_b
         */
        virtual void quadratizeCost2Go(const State&xHat, const Control&uHat, const int&t,
                                       StateMatrix&D,
                                       ControlMatrix&E,
                                       ControlStateMatrix&C,
                                       State&d,
                                       Control&e,
                                       double&scalar_b)
        {
            // ====================== Linearize forward movement model ==================================

            StateMatrix A          = zeros<mat>(xDim, xDim);
            StateControlMatrix B   = zeros<mat>(xDim, uDim);
            State c                = zeros<vec>(xDim);

            const State xHatPrime   = this->system->inverse_dynamics(xHat, uHat);

            linearizeBackward(xHat, uHat, A, B, c);

            // ====================== Quadratize instant cost  ==================================

            ControlStateMatrix P    = zeros<mat>(uDim, xDim);
            StateMatrix Q           = zeros<mat>(xDim, xDim);
            ControlMatrix R         = zeros<mat>(uDim, uDim);

            State q                 = zeros<vec>(xDim);
            Control r               = zeros<vec>(uDim);

            quadratizeInstantCost(xHatPrime, uHat, t, Q, R, P, q, r);

            C   = (B.t() * S.at(t+1) * A) + P;
            D   = (A.t() * S.at(t+1) * A) + Q;
            E   = (B.t() * S.at(t+1) * B) + R;

            d   = A.t() * (s.at(t+1) + S.at(t+1)*c) + q;
            e   = B.t() * (s.at(t+1) + S.at(t+1)*c) + r;

            scalar_b = 0.0; //Not neccesary in SELQR
        }

        /**
         * @brief linearizeForward
         * @param xHat
         * @param uHat
         * @param ABar
         * @param BBar
         * @param cBar
         */
        virtual void linearizeForward(const State&xHat, const Control&uHat,
                                      StateMatrix&ABar, StateControlMatrix&BBar,
                                      State&cBar)
        {
            const State xHatPrime = this->system->move(xHat, uHat);

            ABar    = this->system->jacobianStateInvDynamics(xHatPrime, uHat);
            BBar    = this->system->jacobianControlInvDynamics(xHatPrime, uHat);

            std::cout<<xHat.t()<<std::endl;
            std::cout<<uHat.t()<<std::endl;
            std::cout<<ABar<<std::endl;

            cBar    = xHat - ABar*xHatPrime - BBar*uHat;
        }


        /**
         * @brief linearizeBackward
         * @param xHat
         * @param uHat
         * @param A
         * @param B
         * @param c
         */
        virtual void linearizeBackward(const State&xHat, const Control&uHat,
                                       StateMatrix&A,
                                       StateControlMatrix&B,
                                       State&c)
        {
            const State xHatPrime   = this->system->inverse_dynamics(xHat, uHat);

            A   = this->system->jacobianState(xHatPrime, uHat);
            B   = this->system->jacobianControl(xHatPrime, uHat);
            c   = xHat - A*xHatPrime - B*uHat;

        }


        /**
         * @brief quadratizeInstantCost
         * @param xHat
         * @param uHat
         * @param t
         * @param Q
         * @param R
         * @param P
         * @param q
         * @param r
         */
        virtual void quadratizeInstantCost(const State&xHat, const Control&uHat, const int&t,
                                           StateMatrix&Q,
                                           ControlMatrix&R,
                                           ControlStateMatrix&P,
                                           State&q,
                                           Control&r)
        {

            // Reset all variables
            Q = zeros<mat>(xDim, xDim);
            q = zeros<vec>(xDim);

            R = zeros<mat>(uDim, uDim);
            r = zeros<vec>(uDim);

            P = zeros<mat>(uDim, xDim);

            this->systemCost->quadratize(xHat, uHat, Q, R, P, q, r);

            double unused = 0.0;

            if (t == 0)
            {
                this->initCost->quadratize(xHat, Q, q, unused);
            }

        }

        /**
         * @brief getForwardPolicy
         * @param DBar
         * @param EBar
         * @param CBar
         * @param dBar
         * @param eBar
         * @param scalar
         * @param t
         */
        virtual void getForwardPolicy(const StateMatrix&DBar,
                                      const ControlMatrix&EBar,
                                      const ControlStateMatrix&CBar,
                                      const State&dBar,
                                      const Control&eBar,
                                      const double&scalar,
                                      const int t)
        {

            this->L.at(t) = -solve(EBar, CBar);
            this->l.at(t) = -solve(EBar, eBar);

            const int update = t+1;

            SBar.at(update) = DBar + CBar.t() * this->L.at(t);
            sBar.at(update) = dBar + CBar.t() * this->l.at(t);

            scalar_sBar.at(update) = 0.0;

        }

        /**
         * @brief getBackwardPolicy
         * @param D
         * @param E
         * @param C
         * @param d
         * @param e
         * @param scalar
         * @param t
         */
        virtual void getBackwardPolicy(const StateMatrix&D,
                                       const ControlMatrix&E,
                                       const ControlStateMatrix&C,
                                       const State&d,
                                       const Control&e,
                                       const double&scalar,
                                       const int t)
        {

            this->L.at(t) = -solve(E, C);
            this->l.at(t) = -solve(E, e);

            S.at(t) = D + C.t() * this->L.at(t);
            s.at(t) = d + C.t() * this->l.at(t);

            scalar_s.at(t) = 0.0;
        }


        /**
         * @brief estimateNominalState
         * @param t
         * @param xHat
         */
        void estimateNominalState(const int&t, State&xHat)
        {
            StateMatrix SS    = S.at(t) + SBar.at(t);
            State ss          = s.at(t) + sBar.at(t);

            regularize<xDim>(SS, 1.0e-3, 0.1);
            xHat = - solve(SS, ss);

        }

        /**
         * @brief setInitialConditions
         */
        virtual void setInitialConditions(const std::vector<Control>&nominalU)
        {

        }


        /**
         * @brief reset
         * @param nominalU
         */
        void reset(const std::vector<Control>&nominalU)
        {
            std::fill(S.begin(),        S.end(),        zeros<mat>(xDim, xDim));
            std::fill(s.begin(),        s.end(),        zeros<vec>(xDim));
            std::fill(scalar_s.begin(), scalar_s.end(), 0.0);

            std::fill(SBar.begin(),        SBar.end(),        zeros<mat>(xDim, xDim));
            std::fill(sBar.begin(),        sBar.end(),        zeros<vec>(xDim));
            std::fill(scalar_sBar.begin(), scalar_sBar.end(), 0.0);

            this->clearAll();
            this->setl(nominalU);

        }

        virtual void setParameters()
        {

        }


    protected:

        std::vector< StateMatrix >SBar;
        std::vector< State >sBar;
        std::vector<double>scalar_sBar;


        std::vector< StateMatrix >S;
        std::vector< State>s;
        std::vector<double>scalar_s;


};

#endif // SELQR_HPP
