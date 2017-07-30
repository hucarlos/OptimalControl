#ifndef DDR_HPP
#define DDR_HPP

#include <iostream>
#include <System/BasicSystem.hpp>

class DDR: public BasicSystem<3,2>
{

    public:
        typedef typename BasicSystem<3,2>::State State;
        typedef typename BasicSystem<3,2>::Control Control;

    public:

        DDR(const double&delta_t): dt(delta_t)
        {

        }

        /**
         * @brief Continuous-time dynamics \dot{x} = f(x,u)
         * @param x
         * @param u
         * @return
         */
        inline State continuos_move(const State&x, const Control&u) const
        {
            State xDot = zeros<vec>(3);

            // Differential-drive
            xDot[0] = u[0]*std::cos(x[2]);
            xDot[1] = u[0]*std::sin(x[2]);
            xDot[2] = u[1];


//            xDot[0] = 0.5*(u[0] + u[1])*cos(x[2]);
//            xDot[1] = 0.5*(u[0] + u[1])*sin(x[2]);
//            xDot[2] = (u[1] - u[0])/2.58;

            return xDot;
        }

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

    private:

        const double dt;
};

#endif // DDR_HPP
