#ifndef BASICSYSTEM_HPP
#define BASICSYSTEM_HPP


#include<iostream>
#include <armadillo>

template<uword xDim, uword uDim>
class BasicSystem
{
    public:
        typedef arma::vec::fixed<xDim> State;
        typedef arma::mat::fixed<xDim,xDim> StateMatrix;

        typedef arma::vec::fixed<uDim> Control;

        typedef arma::mat::fixed<xDim,xDim>JacobianState;
        typedef arma::mat::fixed<xDim,uDim>JacobianControl;

    public:
        BasicSystem(const double&delta_t): dt(delta_t), jStep(0.0009765625)
        {

        }

        BasicSystem(const double&delta_t, const double&steapsize): dt(delta_t) ,jStep(steapsize)
        {

        }

        /**
         * @brief operator ()
         * @param state
         * @param control
         * @return
         */
        State operator()(const State&state, const Control&control) const
        {
            return move(state,control);
        }

        /**
         * @brief continuos_move
         * @param x
         * @param u
         * @return
         */
        virtual State continuos_move(const State&x, const Control&u) const = 0;

        /**
         * @brief move Discrete-time dynamics x_{t+1} = g(x_t, u_t)
         * @param x
         * @param u
         * @return
         */
        inline State move(const State& x, const Control& u) const
        {
                State k1  = continuos_move(x, u);
                State k2  = continuos_move(x + 0.5*dt*k1, u);
                State k3  = continuos_move(x + 0.5*dt*k2, u);
                State k4  = continuos_move(x + dt*k3, u);

                return x + (dt/6.0)*(k1 + 2.0*k2 + 2.0*k3 + k4);
        }

        /**
         * @brief Discrete-time inverse dynamics x_t = \bar{g}(x_{t+1}, u_t)
         * @param x
         * @param u
         * @return
         */

        inline State inverse_dynamics(const State& x, const Control& u) const
        {
                State k1 = continuos_move(x, u);
                State k2 = continuos_move(x - 0.5*dt*k1, u);
                State k3 = continuos_move(x - 0.5*dt*k2, u);
                State k4 = continuos_move(x - dt*k3, u);

                return x - (dt/6.0)*(k1 + 2.0*k2 + 2.0*k3 + k4);
        }


        /**
         * @brief Jacobian over state. Numerical if the derivated class has not one
         * @param state
         * @param control
         * @return
         */
        virtual JacobianState jacobianState(const State&state, const Control&control) const
        {
            JacobianState A;
            State ar(state);
            State al(state);

            for (size_t i = 0; i < xDim; ++i)
            {
                ar(i) += jStep;
                al(i) -= jStep;

                A.col(i) = (move(ar, control) - move(al, control)) / (2*jStep);
                ar(i) = al(i) = state(i);
            }

            return A;
        }

        /**
         * @brief Jacobian over control. Numerical if the derivated class has not one
         * @param state
         * @param control
         * @return
         */
        virtual JacobianControl jacobianControl(const State&state, const Control&control) const
        {
            JacobianControl B;
            Control br(control);
            Control bl(control);
            for (size_t i = 0; i < uDim; ++i)
            {
                br(i) += jStep;
                bl(i) -= jStep;
                B.col(i) = (move(state, br) - move(state, bl)) / (2*jStep);
                br(i) = bl(i) = control(i);
            }

            return B;
        }

        /**
         * @brief jacobianStateInvDynamics
         * @param state
         * @param control
         * @return
         */
        virtual JacobianState jacobianStateInvDynamics(const State&state, const Control&control) const
        {
            JacobianState A;
            State ar(state);
            State al(state);

            for (size_t i = 0; i < xDim; ++i)
            {
                ar(i) += jStep;
                al(i) -= jStep;

                A.col(i) = (inverse_dynamics(ar, control) - inverse_dynamics(al, control)) / (2*jStep);
                ar(i) = al(i) = state(i);
            }

            return A;
        }

        /**
         * @brief jacobianControlInvDynamics
         * @param state
         * @param control
         * @return
         */
        virtual JacobianControl jacobianControlInvDynamics(const State&state, const Control&control) const
        {
            JacobianControl B;
            Control br(control);
            Control bl(control);
            for (size_t i = 0; i < uDim; ++i)
            {
                br(i) += jStep;
                bl(i) -= jStep;
                B.col(i) = (inverse_dynamics(state, br) - inverse_dynamics(state, bl)) / (2*jStep);
                br(i) = bl(i) = control(i);
            }

            return B;
        }

    private:
        const double dt;
        const double jStep;
};

#endif // BASICSYSTEM_HPP
