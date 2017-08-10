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

        DDR(const double&delta_t): BasicSystem<3,2>(delta_t)
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

            return xDot;
        }

    private:


};

#endif // DDR_HPP
