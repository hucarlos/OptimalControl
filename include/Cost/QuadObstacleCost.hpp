#ifndef QUADOBSTACLECOST_HPP
#define QUADOBSTACLECOST_HPP

#include <Cost/ObstacleCost.hpp>

template<uword xDim, uword oDim>
class QuadObstacleCost : public ObstaclesCost<xDim, oDim>
{
    public:
        typedef typename Functor<xDim>::Input State;
        typedef typename Functor<xDim>::ImputsMatrix StateMatrix;

    public:

        QuadObstacleCost(): ObstaclesCost<xDim, oDim>()

        {

        }

        /**
         * @brief ObstaclesCost
         * @param robot_radius
         * @param in_obstacles
         */
        QuadObstacleCost(const double&robot_radius, const std::vector< Obstacle<oDim> >&in_obstacles ):
            ObstaclesCost<xDim, oDim>(robot_radius, in_obstacles)
        {
            this->scaleFactor     = 1.0;
            this->obstacleFactor  = 1.0;
        }

        /**
         * @brief ObstaclesCost
         * @param other
         */
        QuadObstacleCost(const ObstaclesCost<xDim, oDim>&other)
        {
            this->copy(other);
        }

        /**
         * @brief evaluate
         * @param state
         * @param control
         * @return
         */
        double evaluate(const State&state) const
        {
            double cost = 0;

            for (size_t i = 0; i < this->obstacles.size(); ++i)
            {
                arma::vec::fixed<oDim>d   =  state.subvec(span(0, oDim-1)) - this->obstacles.at(i).pos;

                d[this->obstacles.at(i).dim] = 0;


                double distr    = std::sqrt(dot(d, d));
                double dist     = distr - this->robot_radius_ - this->obstacles.at(i).radius;
                cost            += this->obstacleFactor * std::exp(-this->scaleFactor*dist);
            }

            for (size_t i = 0; i < oDim; ++i)
            {
                double dist = (state(i) - this->bottomLeft[i]) - this->robot_radius_;
                cost += this->obstacleFactor * std::exp(-this->scaleFactor*dist);
            }

            for (size_t i = 0; i < oDim; ++i)
            {
                double dist = (this->topRight[i] - state(i)) - this->robot_radius_;
                cost += this->obstacleFactor * std::exp(-this->scaleFactor*dist);
            }
            return cost;
        }

        /**
         * @brief quadratizeCost
         * @param x
         * @param Q
         * @param q
         */
        inline void quadratize(const State&x,
                               StateMatrix& Q,
                               State& q,
                               double&scalar) const
        {

            arma::mat::fixed<oDim, oDim> QObs; QObs.zeros();
            arma::vec::fixed<oDim>qObs; qObs.zeros();

            for (size_t i = 0; i < this->obstacles.size(); ++i)
            {
                arma::vec::fixed<oDim> d = x.subvec(0, oDim-1) - this->obstacles.at(i).pos;


                d[this->obstacles.at(i).dim] = 0;

                double distr = std::sqrt(dot(d,d));
                d /= distr;
                double dist = distr - this->robot_radius_ - this->obstacles.at(i).radius;

                vec::fixed<oDim> n          = zeros<vec>(oDim);
                n[this->obstacles[i].dim]   = 1.0;
                vec::fixed<oDim> d_ortho    = skewSymmetric(n)*d;

                double a0 = this->obstacleFactor * std::exp(-this->scaleFactor*dist);
                double a1 = -this->scaleFactor*a0;
                double a2 = -this->scaleFactor*a1;

                double b2 = a1 / distr;

                QObs += a2*(d * d.t()) + b2*(d_ortho * d_ortho.t());
                qObs += a1*d;
            }

            for (size_t i = 0; i < oDim; ++i)
            {
                double dist = (x[i] - this->bottomLeft[i]) - this->robot_radius_;

                arma::vec::fixed<oDim> d;
                d.zeros();

                d[i] = 1.0;

                double a0 = this->obstacleFactor * std::exp(-this->scaleFactor*dist);
                double a1 = -this->scaleFactor*a0;
                double a2 = -this->scaleFactor*a1;

                QObs += a2*(d * d.t());
                qObs += a1*d;
            }

            for (size_t i = 0; i < oDim; ++i)
            {
                double dist = (this->topRight[i] - x[i]) - this->robot_radius_;

                arma::vec::fixed<oDim> d;
                d.zeros();
                d[i] = -1.0;

                double a0 = this->obstacleFactor * exp(-this->scaleFactor*dist);
                double a1 = -this->scaleFactor*a0;
                double a2 = -this->scaleFactor*a1;

                QObs += a2*(d * d.t());
                qObs += a1*d;
            }

            regularize0<oDim>(QObs, 0.0);
            Q.submat(0,0,oDim-1, oDim-1)    += QObs;
            q.subvec(0,oDim-1)              += qObs - QObs * x.subvec(0,oDim-1);

//            std::cout<<QObs<<std::endl;
//            std::cout<<qObs.t()<<std::endl;
            scalar = 0.0;
        }

    private:


};

#endif // QUADOBSTACLECOST_HPP
