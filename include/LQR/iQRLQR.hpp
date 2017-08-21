#ifndef IQRLQR_HPP
#define IQRLQR_HPP

#include <armadillo>
#include <LQR/iLQR.hpp>
#include <Regression/QuadraticRegression.hpp>

template<uword xDim, uword uDim>
class iQRLQR : public iLQR<xDim, uDim>
{
    public:
        typedef typename iLQR<xDim, uDim>::State State;
        typedef typename iLQR<xDim, uDim>::Control Control;

        typedef typename iLQR<xDim, uDim>::StateMatrix StateMatrix;
        typedef typename iLQR<xDim, uDim>::ControlMatrix ControlMatrix;

        typedef typename iLQR<xDim, uDim>::ExtendedState ExtendedState;
        typedef typename iLQR<xDim, uDim>::ExtendedStateMatrix ExtendedStateMatrix;

        typedef typename iLQR<xDim, uDim>::ControlStateMatrix ControlStateMatrix;
        typedef typename iLQR<xDim, uDim>::StateControlMatrix StateControlMatrix;

        typedef typename QuadraticRegression<xDim + uDim>::SAMPLING_MODE SAMPLING_MODE;

    public:
        iQRLQR(const unsigned int ell,
               const BasicSystem<xDim, uDim>*ptr_system,
               const QuadraticCost<xDim> *ptr_init_cost,
               const BasicSystemCost<xDim, uDim>*ptr_system_cost,
               const QuadraticCost<xDim>*ptr_final_cost,
               const bool VIS = true,
               const std::string&name="iQRLQR"):

            iLQR<xDim, uDim>(ell,
                             ptr_system,
                             ptr_init_cost,
                             ptr_system_cost,
                             ptr_final_cost,
                             VIS,
                             name),
            samplingMode(SAMPLING_MODE::GAUSSIAN_S),
            samplingFactor(1.0)
        {
            f_Cost2Go = std::bind(&iQRLQR::cost2Go, this, std::placeholders::_1);

            decreceFactors = 0.5 * arma::ones<vec>(xDim + uDim);

            minEig  = 0.0;
            factEig = 0.1;

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
        void quadratizeCost2Go(const std::vector<State>&xHat,
                               const std::vector<Control>&uHat,
                               StateMatrix&D,
                               ControlMatrix&E,
                               ControlStateMatrix&C,
                               State&d,
                               Control&e,
                               int t)
        {

            ExtendedState mean;
            mean.subvec(0, xDim-1)             = xHat.at(t);
            mean.subvec(xDim, xDim + uDim -1)  = uHat.at(t);

            // ===================== COST TO GO ==================
            setRegression();


            // Results matrix
            ExtendedStateMatrix M;
            ExtendedState m;
            double scalar;

            const double lambda = 0.0;
            radius = initRadius;
            do
            {

                error   = regression.getRegression(f_Cost2Go, mean, M, m, scalar, radius, lambda);
                decreceRadius(radius);
            }
            while(std::abs(error) > epsilon);

            checkM(M);

            m = m - M*mean;

            //Quadratic terms
            C   = M.submat(xDim, 0, size(uDim, xDim));
            D   = M.submat(0,0, size(xDim, xDim));
            E   = M.submat(xDim, xDim, size(uDim, uDim));

            // Linear terms
            d   = m.subvec(0, xDim-1);
            e   = m.subvec(xDim, xDim + uDim - 1);
        }



        /**
         * @brief cost2Go
         * @param state
         * @return
         */
        double cost2Go(const ExtendedState&state) const
        {

            const State x   = state.subvec(0, xDim - 1);
            const Control u = state.subvec(xDim, xDim + uDim - 1);

            const double instantaneous_cost = this->systemCost->evaluate(x, u);

            const State next_state = this->system->move(x, u);

            const double cost_to_go_next = 0.5*dot(next_state, this->S * next_state) + dot(next_state, this->s) + this->scalar_s;

            return instantaneous_cost + cost_to_go_next;
        }

        /**
         * @brief setRegression
         */
        void setRegression()
        {
            regression.setSamplingMode(samplingMode);
            regression.setSamplingFactor(samplingFactor);
            regression.setSeed(seed);
        }

        /**
         * @brief checkM
         * @param M
         */
        void checkM(ExtendedStateMatrix&M)
        {

            toZero<xDim + uDim, xDim + uDim>(M, 1.0e-9);
            regularize<xDim + uDim>(M, minEig, factEig);
        }

        /**
         * @brief decreceRadius
         * @param radius
         */
        void decreceRadius(ExtendedState&radius)
        {

//            for(unsigned int  r=0; r<(xDim + uDim); r++)
//            {
//                if(r == 2)
//                {
//                    radius(r) *= 0.1;
//                }
//                else
//                {
//                    radius(r) *= 0.5;
//                }
//            }

            radius.subvec(0, 2)     *= 0.5;
            radius.subvec(3, 5)     *= 0.5;
            radius.subvec(6, 8)     *= 0.5;
            radius.subvec(9, 11)    *= 0.5;
            radius.subvec(12, 15)   *= 0.5;

            radius %= decreceFactors;
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
         * @brief getRadius
         * @return
         */
        ExtendedState getRadius() const
        {
            return radius;
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
        void setEpsilon(double value)
        {
            epsilon = value;
        }

        /**
         * @brief getGaussian_sampling
         * @return
         */
        SAMPLING_MODE getSamplingMode() const
        {
            return samplingMode;
        }

        /**
         * @brief setGaussian_sampling
         * @param value
         */
        void setSamplingMode(SAMPLING_MODE value)
        {
            samplingMode = value;
        }

        /**
         * @brief getSampling_factor
         * @return
         */
        unsigned int getSamplingFactor() const
        {
            return samplingFactor;
        }

        /**
         * @brief setSampling_factor
         * @param value
         */
        void setSamplingFactor(unsigned int value)
        {
            samplingFactor = value;
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

    protected:

        std::function< double(const ExtendedState&)> f_Cost2Go;

        // Regression quadratic parameters
        ExtendedState initRadius;
        ExtendedState radius;
        double epsilon;
        double error;

        arma::vec::fixed<xDim+uDim>decreceFactors;

        QuadraticRegression<xDim + uDim>regression;
        SAMPLING_MODE samplingMode;
        unsigned int samplingFactor;
        unsigned int seed;

        double minEig;
        double factEig;


};


#endif // IQRLQR_HPP
