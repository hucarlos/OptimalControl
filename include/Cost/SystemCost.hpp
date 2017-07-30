#ifndef SYSTEMCOST_HPP
#define SYSTEMCOST_HPP

#include <iostream>

#include <Cost/BasicSystemCost.hpp>
#include <Cost/QuadraticCost.hpp>
#include <Cost/ObstacleCost.hpp>

template<uword xDim, uword uDim>
class SystemCost : public BasicSystemCost<xDim, uDim>
{
    public:

        typedef typename BasicSystemCost<xDim, uDim>::State State;
        typedef typename BasicSystemCost<xDim, uDim>::Control Control;

        typedef typename BasicSystemCost<xDim, uDim>::StateMatrix StateMatrix;
        typedef typename BasicSystemCost<xDim, uDim>::ControlMatrix ControlMatrix;

        typedef typename BasicSystemCost<xDim, uDim>::ControlStateMatrix ControlStateMatrix;

    public:
        SystemCost(): BasicSystemCost<xDim, uDim>(),
            control_cost(NULL),
            obstacles_cost(NULL)

        {

        }

        ~SystemCost()
        {
            control_cost    = NULL;
            obstacles_cost  = NULL;
        }

        SystemCost(QuadraticCost<uDim>*ptr_control_cost,
                   ObstaclesCost<xDim,2>*ptr_obstacles_cost):

              BasicSystemCost<xDim, uDim>(),
              control_cost(ptr_control_cost),
              obstacles_cost(ptr_obstacles_cost)

        {

        }

        double evaluate(const State&state, const Control&control) const
        {

            double cost = 0.0;

            // Control cost
            cost += control_cost->evaluate(control);

            // Obstacle cost
            cost += obstacles_cost->evaluate(state);


            return cost;
        }

        /**
         * @brief operator ()
         * @param state
         * @return
         */
        double evaluateState(const State&state) const
        {
            return this->obstacles_cost->evaluate(state);
        }


        /**
         * @brief quadratize
         * @param state
         * @param control
         * @param Qt
         * @param Rt
         * @param Pt
         * @param qt
         * @param rt
         */
        void quadratize(const State&state, const Control&control,
                        StateMatrix&Qt, ControlMatrix&Rt, ControlStateMatrix&Pt,
                        State&qt, Control&rt) const
        {

            double scalarC      = 0.0;
            double scalarObs    = 0.0;

            control_cost->quadratize(control, Rt, rt, scalarC);

            obstacles_cost->quadratize(state, Qt, qt, scalarObs);

            Pt = zeros<mat>(uDim, xDim);

        }

    protected:

        const QuadraticCost<uDim> *control_cost;
        const ObstaclesCost<xDim,2> *obstacles_cost;
};



#endif // SYSTEMCOST_HPP
