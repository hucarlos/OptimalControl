#ifndef OBSTACLECOST_HPP
#define OBSTACLECOST_HPP

#include <iostream>
#include <vector>

#include<Cost/Functor.hpp>
#include <Cost/Obstacle.hpp>

#include <Utils/Utils.hpp>

/**
 * @brief The ObstaclesCost class
 */
template<size_t xDim, size_t obsDim>
class ObstaclesCost: public Functor<xDim>
{
    public:
        typedef typename Functor<xDim>::Input State;
        typedef typename Functor<xDim>::ImputsMatrix StateMatrix;

    public:

        ObstaclesCost():
            Functor<xDim>(),
            robot_radius_(0.0),
            scaleFactor(1.0),
            obstacleFactor(1.0)
        {

        }

        /**
         * @brief ObstaclesCost
         * @param robot_radius
         * @param in_obstacles
         */
        ObstaclesCost(const double&robot_radius, const std::vector< Obstacle<obsDim> >&in_obstacles ):
            Functor<xDim>(), obstacles(in_obstacles), robot_radius_(robot_radius)
        {
            scaleFactor     = 1.0;
            obstacleFactor  = 1.0;
        }

        /**
         * @brief ObstaclesCost
         * @param other
         */
        ObstaclesCost(const ObstaclesCost<xDim, obsDim>&other):
            Functor<xDim>()
        {
            this->copy(other);
        }

        void copy(const ObstaclesCost<xDim, obsDim>&other)
        {
            robot_radius_   = other.robot_radius_;
            scaleFactor     = other.scaleFactor;
            obstacleFactor  = other.obstacleFactor;

            obstacles       = other.obstacles;

            bottomLeft      = other.bottomLeft;
            topRight        = other.topRight;
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

            for (size_t i = 0; i < obstacles.size(); ++i)
            {
                arma::vec::fixed<obsDim>d   =  state.subvec(span(0, obsDim-1)) - obstacles.at(i).pos;

                double distr    = std::sqrt(dot(d, d));
                double dist     = distr - robot_radius_ - obstacles.at(i).radius;
                cost            += obstacleFactor * std::exp(-scaleFactor*dist);
            }

            for (size_t i = 0; i < obsDim; ++i)
            {
                double dist = (state(i) - bottomLeft[i]) - robot_radius_;
                cost += obstacleFactor * std::exp(-scaleFactor*dist);
            }

            for (size_t i = 0; i < obsDim; ++i)
            {
                double dist = (topRight[i] - state(i)) - robot_radius_;
                cost += obstacleFactor * std::exp(-scaleFactor*dist);
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

            arma::mat::fixed<obsDim, obsDim> QObs; QObs.zeros();
            arma::vec::fixed<obsDim>qObs; qObs.zeros();

            for (size_t i = 0; i < obstacles.size(); ++i)
            {
                arma::vec::fixed<obsDim> d = x.subvec(0, obsDim-1) - obstacles.at(i).pos;

                double distr = std::sqrt(dot(d,d));
                d /= distr;
                double dist = distr - robot_radius_ - obstacles.at(i).radius;

                arma::vec::fixed<obsDim> d_ortho;
                d_ortho[0] = d[1];
                d_ortho[1] = -d[0];

                double a0 = obstacleFactor * std::exp(-scaleFactor*dist);
                double a1 = -scaleFactor*a0;
                double a2 = -scaleFactor*a1;

                double b2 = a1 / distr;

                QObs += a2*(d * d.t()) + b2*(d_ortho * d_ortho.t());
                qObs += a1*d;
            }

            for (size_t i = 0; i < obsDim; ++i)
            {
                double dist = (x[i] - bottomLeft[i]) - robot_radius_;

                arma::vec::fixed<obsDim> d;
                d.zeros();

                d[i] = 1.0;

                double a0 = obstacleFactor * std::exp(-scaleFactor*dist);
                double a1 = -scaleFactor*a0;
                double a2 = -scaleFactor*a1;

                QObs += a2*(d * d.t());
                qObs += a1*d;
            }

            for (size_t i = 0; i < obsDim; ++i)
            {
                double dist = (topRight[i] - x[i]) - robot_radius_;

                arma::vec::fixed<obsDim> d;
                d.zeros();
                d[i] = -1.0;

                double a0 = obstacleFactor * exp(-scaleFactor*dist);
                double a1 = -scaleFactor*a0;
                double a2 = -scaleFactor*a1;

                QObs += a2*(d * d.t());
                qObs += a1*d;
            }

            regularize0<obsDim>(QObs, 0.1);
            Q.submat(0,0,obsDim-1, obsDim-1)    += QObs;
            q.subvec(0,obsDim-1)                += qObs - QObs * x.subvec(0,obsDim-1);

            scalar = 0.0;
        }

        // ====================================== properties getter and setter===================================

        double getRobotRadius() const
        {
            return robot_radius_;
        }

        void setRobotRadius(const double&robot_radius)
        {
            robot_radius_ = robot_radius;
        }

        std::vector< Obstacle<xDim> > getObstacles() const
        {
            return obstacles;
        }

        arma::vec::fixed<obsDim> getBottomLeft() const
        {
            return bottomLeft;
        }

        void setBottomLeft(const arma::vec::fixed<obsDim> &value)
        {
            bottomLeft = value;
        }

        arma::vec::fixed<obsDim> getTopRight() const
        {
            return topRight;
        }

        void setTopRight(const arma::vec::fixed<obsDim> &value)
        {
            topRight = value;
        }

        double getScaleFactor() const
        {
            return scaleFactor;
        }

        void setScaleFactor(double value)
        {
            scaleFactor = value;
        }

        double getObstacleFactor() const
        {
            return obstacleFactor;
        }

        void setObstacleFactor(double value)
        {
            obstacleFactor = value;
        }

    protected:

        std::vector< Obstacle<obsDim> > obstacles;

        double robot_radius_;

        arma::vec::fixed<obsDim> bottomLeft;
        arma::vec::fixed<obsDim> topRight;

        double scaleFactor;
        double obstacleFactor;

};

#endif // OBSTACLECOST_HPP
