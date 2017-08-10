#ifndef QUADROTOR_HPP
#define QUADROTOR_HPP

#include <System/BasicSystem.hpp>
#include <Utils/Utils.hpp>

class Quadrotor: public BasicSystem<12,4>
{

    public:
        typedef typename BasicSystem<12,4>::State State;
        typedef typename BasicSystem<12,4>::Control Control;

    public:

        Quadrotor(const double&delta_t): BasicSystem<12,4>(delta_t)
        {
            gravity = 9.80665;                  // gravity,   m/s^2
            eX(0) = 1; eX(1) = 0; eX(2) = 0;
            eY(0) = 0; eY(1) = 1; eY(2) = 0;
            eZ(0) = 0; eZ(1) = 0; eZ(2) = 1;

            // Quadrotor parameters
            mass        = 0.5;                  // mass, kg  (source: paper)
            inertia     = 0.05 * eye<mat>(3,3);   // moment of inertia matrix
            momentConst = 1.5e-9 / 6.11e-8;     // ratio between force and moment of rotor
            dragConst   = 0.15;                 //0.15;
            length      = 0.3429/2;             // distance between center and rotor, m
            invInertia  = !inertia;
        }

        /**
         * @brief Continuous-time dynamics \dot{x} = f(x,u)
         * @param x
         * @param u
         * @return
         */
        inline State continuos_move(const State&x, const Control&u) const
        {
            State xDot = zeros<vec>(12);

            vec::fixed<3> p = x.subvec(0, 2);
            vec::fixed<3> v = x.subvec(3, 5);
            vec::fixed<3> r = x.subvec(6, 8);
            vec::fixed<3> w = x.subvec(9, 11);

            // \dot{p} = v
            xDot.subvec(0,2) = v;

            // \dot{v} = [0,0,-g]^T + R*exp([r])*[0,0,(f_1 + f_2 + f_3 + f_4) / m]^T;
            vec::fixed<3>temp = -gravity*eZ + exp(skewSymmetric(r)) * ((u[0]+u[1]+u[2]+u[3])/mass)*eZ - dragConst*v/mass;
            xDot.subvec(3, 5) = temp;

            // \dot{r} = w + 0.5*skewSymmetric(r)*w + (1.0/tr(~r*r))*(1.0 - 0.5*sqrt(tr(~r*r))/tan(0.5*sqrt(tr(~r*r))))*skewSymmetric(r)*(skewSymmetric(r)*w)
            double l = std::sqrt(r[0]*r[0]+r[1]*r[1]+r[2]*r[2]);
            if (0.5*l > 0.0)
            {
                xDot.subvec(6, 8) =  w + 0.5*skewSymmetric(r)*w
                                       + (1.0 - 0.5*l/std::tan(0.5*l))*skewSymmetric(r / l)*(skewSymmetric(r / l)*w);
            }
            else
            {
                xDot.subvec(6, 8) =  w;
            }

            // \dot{w} = J^{-1}*([l*(f_2 - f_4), l*(f_3 - f_1), (f_1 - f_2 + f_3 - f_4)*k_M]^T - [w]*J*w)
            xDot.subvec(9, 11) = invInertia*( length*(u[1] - u[3])*eX + length*(u[2] - u[0])*eY
                                 + (u[0] - u[1] + u[2] - u[3])*momentConst*eZ
                                 - skewSymmetric(w)*inertia*w);

            return xDot;
        }

        /**
         * @brief getGravity
         * @return
         */
        double getGravity() const
        {
            return gravity;
        }

        /**
         * @brief getMass
         * @return
         */
        double getMass() const
        {
            return mass;
        }

    protected:

        double gravity;         // gravity,   m/s^2
        vec::fixed<3> eX, eY, eZ;      // unit vectors

        // quadrotor constants
        double mass;                // mass, kg
        mat::fixed<3,3> inertia;    // moment of inertia matrix
        mat::fixed<3,3> invInertia; // moment of inertia matrix
        double momentConst;         // ratio between force and moment of rotor
        double dragConst;           // ratio between speed and opposite drag force
        double length;              // distance between center and rotor, m
        double minForce;            // minimum force of rotor, N
        double maxForce;            // maximum force of rotor, N
};

#endif // QUADROTOR_HPP
