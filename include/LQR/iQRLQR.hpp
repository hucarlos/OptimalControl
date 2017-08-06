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
                             name)
        {
            f_Cost2Go = std::bind(&iQRLQR::cost2Go, this, std::placeholders::_1);

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

            QuadraticRegression<xDim + uDim>regression;

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

            toZero<xDim + uDim, xDim + uDim>(M, 1.0e-9);
            const double max = M.diag().max();//0.1;
            regularize<xDim + uDim>(M, 0.0, max);
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
         * @brief decreceRadius
         * @param radius
         */
        void decreceRadius(ExtendedState&radius)
        {

            for(int r=0; r<(xDim + uDim); r++)
            {
                if(r == 2)
                {
                    radius(r) *= 0.1;
                }
                else
                {
                    radius(r) *= 0.5;
                }
            }
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

    protected:

        std::function< double(const ExtendedState&)> f_Cost2Go;

        // Regression quadratic parameters
        ExtendedState initRadius;
        ExtendedState radius;
        double epsilon;
        double error;


};


#endif // IQRLQR_HPP
