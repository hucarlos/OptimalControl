#ifndef IQRSELQR_HPP
#define IQRSELQR_HPP


#include <LQR/SELQR.hpp>
#include <Regression/QuadraticRegression.hpp>

template<uword xDim, uword uDim>
class iQRSELQR : public SELQR<xDim, uDim>
{
    public:
        typedef typename SELQR<xDim,uDim>::State State;
        typedef typename SELQR<xDim,uDim>::Control Control;

        typedef typename SELQR<xDim,uDim>::StateMatrix StateMatrix;
        typedef typename SELQR<xDim,uDim>::ControlMatrix ControlMatrix;

        typedef typename SELQR<xDim,uDim>::ExtendedState ExtendedState;
        typedef typename SELQR<xDim,uDim>::ExtendedStateMatrix ExtendedStateMatrix;

        typedef typename SELQR<xDim,uDim>::ControlStateMatrix ControlStateMatrix;
        typedef typename SELQR<xDim,uDim>::StateControlMatrix StateControlMatrix;

        typedef typename QuadraticRegression<xDim + uDim>::SAMPLING_MODE SAMPLING_MODE;

        iQRSELQR(const unsigned int ell,
                 const BasicSystem<xDim, uDim>*ptr_system,
                 const QuadraticCost<xDim> *ptr_init_cost,
                 const BasicSystemCost<xDim, uDim>*ptr_system_cost,
                 const QuadraticCost<xDim>*ptr_final_cost,
                 const bool VIS = true,
                 const std::string&name="iQRSELQR"):

            SELQR<xDim, uDim>(ell, ptr_system, ptr_init_cost, ptr_system_cost, ptr_final_cost, VIS, name),
            epsilon(1.0e-2)
        {
            decreceFactors = 0.5 * arma::ones<vec>(xDim + uDim);

            minEig  = 0.0;
            factEig = 0.1;
        }

        void quadratizeCost2Come(const State&xHat, const Control&uHat, const int&t,
                                 StateMatrix&D,
                                 ControlMatrix&E,
                                 ControlStateMatrix&C,
                                 State&d,
                                 Control&e,
                                 double&scalar_b)
        {

            const State xHatPrime   = this->system->move(xHat, uHat);

            ExtendedState mean;
            mean.subvec(0, xDim-1)                 = xHatPrime;
            mean.subvec(xDim, xDim + uDim - 1)     = uHat;

            // Results matrix
            ExtendedStateMatrix M;
            ExtendedState m;
            double scalar;

            std::function< double(const ExtendedState&)> function;

            if(t==0)
            {
                function = std::function< double(const ExtendedState&)>(std::bind(&iQRSELQR<xDim, uDim>::cost2Come0,
                                                                                  this, std::placeholders::_1));
            }
            else
            {
                function = std::function< double(const ExtendedState&)>(std::bind(&iQRSELQR<xDim, uDim>::cost2Come,
                                                                                  this, std::placeholders::_1));
            }

            const double evalEpsilon = estimateEpsilon();

            const double lambda = 0.0;

            radius = initRadius;
            do
            {

                error   = regression.getRegression(function, mean, M, m, scalar, radius, lambda);
                decreceRadius(radius);
            }
            while(std::abs(error) > evalEpsilon);


            checkM(M);

            m = m - M*mean;

            //Quadratic terms
            C   = M.submat(xDim, 0, size(uDim, xDim));
            D   = M.submat(0,0, size(xDim, xDim));
            E   = M.submat(xDim, xDim, size(uDim, uDim));

            // Linear terms
            d   = m.subvec(0, xDim - 1);
            e   = m.subvec(xDim, xDim + uDim -1);

            scalar_b = scalar;
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
        void quadratizeCost2Go(const State&xHat, const Control&uHat, const int&t,
                               StateMatrix&D,
                               ControlMatrix&E,
                               ControlStateMatrix&C,
                               State&d,
                               Control&e,
                               double&scalar_b)
        {

            const State xHatPrime   = this->system->inverse_dynamics(xHat, uHat);

            ExtendedState mean;
            mean.subvec(0, xDim - 1)               = xHatPrime;
            mean.subvec(xDim, xDim + uDim -1)      = uHat;

            // Results matrix
            ExtendedStateMatrix M;
            ExtendedState m;
            double scalar;

            std::function< double(const ExtendedState&)> function;

            if(t==0)
            {
                function = std::function< double(const ExtendedState&)>(std::bind(&iQRSELQR<xDim, uDim>::cost2Go0,
                                                                                  this, std::placeholders::_1));
            }
            else
            {
                function = std::function< double(const ExtendedState&)>(std::bind(&iQRSELQR<xDim, uDim>::cost2Go,
                                                                                  this, std::placeholders::_1));
            }

            const double evalEpsilon = estimateEpsilon();

            const double lambda = 0.0;
            radius = initRadius;
            do
            {

                error   = regression.getRegression(function, mean, M, m, scalar, radius, lambda);
                decreceRadius(radius);
            }
            while(std::abs(error) > evalEpsilon);


            checkM(M);

            m = m - M*mean;


            //Quadratic terms
            C   = M.submat(xDim, 0, size(uDim, xDim));
            D   = M.submat(0,0, size(xDim, xDim));
            E   = M.submat(xDim, xDim, size(uDim, uDim));


            // Linear terms
            d   = m.subvec(0, xDim - 1);
            e   = m.subvec(xDim, xDim + uDim -1);

            scalar_b = scalar;
        }


        /**
         * @brief checkM
         * @param M
         */
        void checkM(ExtendedStateMatrix&M)
        {

            toZero<xDim + uDim, xDim + uDim>(M, 1.0e-9);

            if(this->iters2 < 10)
            {
                regularize<xDim + uDim>(M, minEig, 100);
            }
            else
            {
                regularize<xDim + uDim>(M, minEig, factEig);
            }

        }


        /**
         * @brief setParameters
         */
        void setParameters()
        {

        }


        /**
         * @brief setInitialConditions
         * @param start
         * @param initStart
         * @param newCost
         * @param oldCost
         */
        void setInitialConditions(const std::vector<Control>&nominalU)
        {
            epsilon *= 0.1;
            estimateEpsilon();

            initRadius *= 0.5;
//            initRadius %= decreceFactors;

            if(this->vis)
                std::cout<<"Radius :"<<initRadius.t();
            
            this->reset(nominalU);
        }

        /**
         * @brief cost2Come
         * @param state
         * @return
         */
        double cost2Come(const ExtendedState&state) const
        {
            const unsigned int t = this->time;
            const StateMatrix St = this->SBar.at(t);
            const State st       = this->sBar.at(t);
            const double scalar  = this->scalar_sBar.at(t);

            const State x   = state.subvec(0, xDim-1);
            const Control u = state.subvec(xDim, xDim + uDim -1);

            const State xPrime = this->system->inverse_dynamics(x, u);

            const double instantaneous_cost = this->systemCost->evaluate(xPrime, u);

            const double cost_to_go_backward = 0.5*dot( xPrime, St*xPrime) + dot(xPrime, st) + scalar;

            return instantaneous_cost + cost_to_go_backward;
        }


        /**
         * @brief cost2Go
         * @param state
         * @return
         */
        double cost2Go(const ExtendedState&state) const
        {

            const unsigned int t = this->time + 1;
            const StateMatrix St = this->S.at(t);
            const State st       = this->s.at(t);
            const double scalar  = this->scalar_s.at(t);

            const State x   = state.subvec(0, xDim - 1);
            const Control u = state.subvec(xDim, xDim + uDim -1);

            const double instantaneous_cost = this->systemCost->evaluate(x, u);

            const State next_state = this->system->move(x, u);

            const double cost_to_go_next = 0.5*dot(next_state, St*next_state) + dot(next_state, st) + scalar;

            return instantaneous_cost + cost_to_go_next;
        }

        /**
         * @brief cost2Come0
         * @param state
         * @return
         */
        double cost2Come0(const ExtendedState&state) const
        {

            const unsigned int t = 0;
            const StateMatrix St = this->SBar.at(t);
            const State st       = this->sBar.at(t);
            const double scalar  = this->scalar_sBar.at(t);

            const State x   = state.subvec(0, xDim-1);
            const Control u = state.subvec(xDim, xDim + uDim -1);

            const State xPrime = this->system->inverse_dynamics(x, u);

            const double instantaneous_cost = this->systemCost->evaluate(xPrime, u) + this->initCost->evaluate(xPrime);

            const double cost_to_go_backward = 0.5*dot( xPrime, St*xPrime) + dot(xPrime, st) + scalar;

            return instantaneous_cost + cost_to_go_backward;
        }

        /**
         * @brief cost2Go0
         * @param state
         * @return
         */
        double cost2Go0(const ExtendedState&state) const
        {
            const unsigned int t = 1;
            const StateMatrix St = this->S.at(t);
            const State st       = this->s.at(t);
            const double scalar  = this->scalar_s.at(t);

            const State x   = state.subvec(0, xDim - 1);
            const Control u = state.subvec(xDim, xDim + uDim -1);

            const double instantaneous_cost = this->systemCost->evaluate(x, u) + this->initCost->evaluate(x);

            const State next_state = this->system->move(x, u);

            const double cost_to_go_next = 0.5*dot(next_state, St*next_state) + dot(next_state, st) + scalar;

            return instantaneous_cost + cost_to_go_next;

        }

        /**
         * @brief estimateEpsilon
         * @return
         */
        double estimateEpsilon()
        {

            epsilon = std::max(1.0e-10, epsilon);
            return epsilon;
        }

        /**
         * @brief selectRadius
         * @param state
         * @param radius
         */
        void decreceRadius(ExtendedState&radius)
        {

            radius %= decreceFactors;
            
        }


        /**
         * @brief printConfiguration
         */
        void printConfiguration()
        {
            std::cout   << this->iterations()   <<' '   << this->accum      << ' '
                        << this->epsilon        <<' '
                        << error                <<' '   << initRadius.t()   << ' '
                        <<std::endl;
        }


        /**
         * @brief getInitRadius
         * @return
         */
        ExtendedState getInitRadius() const
        {
            return initRadius;
        }

        /**
         * @brief setInitRadius
         * @param value
         */
        void setInitRadius(const ExtendedState &value)
        {
            initRadius = value;
        }

        /**
         * @brief getEpsilon
         * @return
         */
        double getEpsilon() const
        {
            return epsilon;
        }

        /**
         * @brief setEpsilon
         * @param value
         */
        void setEpsilon(const double&value)
        {
            epsilon = value;
        }

        /**
         * @brief setDecreceFactors
         * @param input
         */
        void setDecreceFactors(arma::vec::fixed<xDim+uDim>&input)
        {
            decreceFactors = input;
        }

        /**
         * @brief setDecreceFactors
         * @param factor
         */
        void setDecreceFactors(const double&factor)
        {
            decreceFactors = factor * arma::ones<vec>(xDim + uDim);
        }

        /**
         * @brief setMinEig
         * @param min
         */
        void setMinEig(const double&min)
        {
            minEig = min;
        }

        /**
         * @brief setFactEig
         * @param fac
         */
        void setFactEig(const double&fac)
        {
            factEig = fac;
        }

        // For regression properties
        /**
         * @brief setGaussian_sampling
         * @param value
         */
        void setSamplingMode(SAMPLING_MODE value)
        {
            regression.setSamplingMode(value);
        }

        /**
         * @brief setSampling_factor
         * @param value
         */
        void setSamplingFactor(unsigned int value)
        {
           regression.setSamplingFactor(value);
        }

        /**
         * @brief setParallel
         * @param parallel
         */
        void setParallel(bool parallel)
        {
            regression.setParallel(parallel);
        }

        /**
         * @brief setRegressionSeed
         * @param seed
         */
        void setRegressionSeed(const unsigned int&seed)
        {
            regression.setSeed(seed);
        }


protected:

        // Regression quadratic parameters
        ExtendedState initRadius;
        ExtendedState radius;
        double epsilon;
        double error;

        QuadraticRegression<xDim + uDim>regression;


        arma::vec::fixed<xDim+uDim>decreceFactors;

        double minEig;
        double factEig;


};

#endif // IQRSELQR_HPP
