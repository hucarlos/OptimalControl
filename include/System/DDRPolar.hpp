#ifndef DDRPOLAR_HPP
#define DDRPOLAR_HPP

#include <iostream>
#include <BasicSystem.hpp>

using namespace Eigen;
using namespace std;

template<size_t xDim>
/**
 * @brief toCartesian
 * @param polar
 * @param cartesian
 */
inline void toCartesian(const Matrix<double, xDim,1>&polar,
                        Matrix<double, xDim,1>&cartesian)
{
    cartesian[0] = polar[0] * std::cos(polar[1]);
    cartesian[1] = polar[0] * std::sin(polar[1]);
    cartesian[2] = polar[1] - polar[2] + M_PI;
    cartesian[3] = polar[3];
}


template<size_t xDim>
inline void toPolar(const Matrix<double, xDim,1>&cartesian,
                    Matrix<double, xDim,1>&polar)
{
    polar[0] = std::sqrt(pow(cartesian[0],2) + pow(cartesian[1],2));
    polar[1] = std::atan2(cartesian[1], cartesian[0]);
    polar[2] = polar[1] - cartesian[2] + M_PI;
    polar[3] = cartesian[3];
}

/**
 * @brief The DDRPolar class
 */
class DDRPolar: public BasicSystem<4,2>
{


    public:

        DDRPolar()
        {

        }

        /**
         * @brief Continuous-time dynamics \dot{x} = f(x,u)
         * @param x
         * @param u
         * @return
         */
        inline BasicSystem<4,2>::State continuos_move(const BasicSystem<4,2>::State&x,
                                                      const BasicSystem<4,2>::Control&u) const
        {
            BasicSystem<4,2>::State xDot = BasicSystem<4,2>::State::Zero();

            // Differential-drive
            const double fac    = 1.0/x[0];
            const double phi    = x[2];
            const double cosine = std::cos(phi);
            const double sine   = std::sin(phi);

            xDot[0] = -cosine    * u[0];
            xDot[1] = fac * sine * u[0];
            xDot[2] = (fac * sine* u[0]) - u[1];

            xDot[3] = 0.0;
            return xDot;
        }

        /**
         * @brief move Discrete-time dynamics x_{t+1} = g(x_t, u_t)
         * @param x
         * @param u
         * @return
         */
        inline BasicSystem<4,2>::State move(const BasicSystem<4,2>::State& x,
                                            const BasicSystem<4,2>::Control& u) const
        {
                const double dt = std::exp(x[3]);

                BasicSystem<4,2>::State k1  = continuos_move(x, u);
                BasicSystem<4,2>::State k2  = continuos_move(x + 0.5*dt*k1, u);
                BasicSystem<4,2>::State k3  = continuos_move(x + 0.5*dt*k2, u);
                BasicSystem<4,2>::State k4  = continuos_move(x + dt*k3, u);

                return x + (dt/6.0)*(k1 + 2.0*k2 + 2.0*k3 + k4);
        }

        /**
         * @brief Discrete-time inverse dynamics x_t = \bar{g}(x_{t+1}, u_t)
         * @param x
         * @param u
         * @return
         */

        inline BasicSystem<4,2>::State inverse_dynamics(const BasicSystem<4,2>::State& x,
                                                        const BasicSystem<4,2>::Control& u) const
        {
                const double dt = std::exp(x[3]);

                BasicSystem<4,2>::State k1 = continuos_move(x, u);
                BasicSystem<4,2>::State k2 = continuos_move(x - 0.5*dt*k1, u);
                BasicSystem<4,2>::State k3 = continuos_move(x - 0.5*dt*k2, u);
                BasicSystem<4,2>::State k4 = continuos_move(x - dt*k3, u);

                return x - (dt/6.0)*(k1 + 2.0*k2 + 2.0*k3 + k4);
        }




    private:
};

#endif // DDRPOLAR_HPP
