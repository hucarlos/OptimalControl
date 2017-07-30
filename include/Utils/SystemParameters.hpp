#ifndef SYSTEMPARAMETERS_HPP
#define SYSTEMPARAMETERS_HPP

#include <iostream>
#include <armadillo>

#include<yaml-cpp/yaml.h>

using namespace std;
using namespace arma;

template<int xDim, int uDim>
class ObstaclesParameters
{
    public:

        typedef vec::fixed<xDim>State;
        typedef vec::fixed<uDim>Control;

        typedef mat::fixed<xDim, xDim>StateMatrix;
        typedef mat::fixed<uDim, uDim>ControlMatrix;

        typedef vec::fixed<xDim + uDim>ExtendedState;

    public:

        State initState;
        State goalState;

        Control nominalControl;

        unsigned int horizont;
        unsigned int maxIters;
        double epsilon;

        double controlFactor;
        double stateFactor;

        double obstacleFactor;
        double scaleFactor;
        double robotRadius;
    
    
        struct lqr
        {
            ExtendedState initRadius;
            double epsilon;
            
        } iQRLQR, iQRSELQR;

        string mapName;

};

template<int xDim, int uDim>
class ObstacleSystemParameters
{

    public:

        /**
         * @brief loadYAMLParams
         * @param filename
         * @param parameters
         */
        void loadYAMLParams(const std::string&filename,  ObstaclesParameters<xDim, uDim>&parameters)
        {
            try
            {

                YAML::Node parentNode = YAML::LoadFile(filename);

                for (auto it_parent = parentNode.begin(); it_parent != parentNode.end(); it_parent++)
                {
                    YAML::Node key      = it_parent->first;
                    YAML::Node child    = it_parent->second;

                    const std::string str_key   = key.as<std::string>();


                    if (parentNode.IsMap())
                    {

                        if(child.IsMap())
                        {

                            // Now that you've got a map, you can iterate through it:
                            for (auto it = child.begin(); it != child.end(); ++it)
                            {


                                const YAML::Node childKey     = it->first;
                                const YAML::Node childValue   = it->second;

                                if(str_key == "xStart")
                                {
                                    parameters.initState[0]    = childValue[0].as<double>();
                                    parameters.initState[1]    = childValue[1].as<double>();
                                    parameters.initState[2]    = (M_PI/180)*childValue[2].as<double>();
                                }

                                else if(str_key == "xGoal")
                                {
                                    parameters.goalState[0]   = childValue[0].as<double>();
                                    parameters.goalState[1]   = childValue[1].as<double>();
                                    parameters.goalState[2]   = (M_PI/180.0)*childValue[2].as<double>();

                                }

                                else if(str_key == "nominalControl")
                                {
                                    parameters.nominalControl[0]    = childValue[0].as<double>();
                                    parameters.nominalControl[1]    = childValue[1].as<double>();

                                }

                                else if(str_key == "quadraticFactors")
                                {
                                    parameters.stateFactor      = childValue[0].as<double>();
                                    parameters.controlFactor    = childValue[1].as<double>();

                                }

                                else if(str_key == "obstaclesFactors")
                                {
                                    parameters.obstacleFactor    = childValue[0].as<double>();
                                    parameters.scaleFactor       = childValue[1].as<double>();
                                    parameters.robotRadius       = childValue[2].as<double>();

                                }

                                else if(str_key == "paramsLQR")
                                {
                                    parameters.horizont     = childValue[0].as<unsigned int>();
                                    parameters.epsilon      = childValue[1].as<double>();
                                    parameters.maxIters     = childValue[2].as<unsigned int>();

                                }

                                else if(str_key == "mapFile")
                                {
                                    parameters.mapName    = childValue[0].as<std::string>();
                                }


                                else if(str_key == "iQRLQR")
                                {

                                    const std::string id    = childKey.as<std::string>();

                                    if(id=="epsilon")
                                    {
                                        parameters.iQRLQR.epsilon   = childValue[0].as<double>();
                                    }

                                    if(id == "radius")
                                    {

                                        for(int i=0; i< (xDim + uDim); i++)
                                        {
                                            parameters.iQRLQR.initRadius(i) = childValue[i].as<double>();
                                        }

                                    }

                                }

                                else if(str_key == "iQRSELQR")
                                {

                                    const std::string id    = childKey.as<std::string>();

                                    if(id=="epsilon")
                                    {
                                        parameters.iQRSELQR.epsilon = childValue[0].as<double>();
                                    }

                                    if(id == "radius")
                                    {
                                        for(int i=0; i< (xDim + uDim); i++)
                                        {
                                           parameters.iQRSELQR.initRadius(i) = childValue[i].as<double>();
                                        }

                                    }                                    
                                }
                                
                                else
                                {
                                    std::string errorMessage = std::string("Error: on read file: ") + str_key + (" not know");
                                    throw std::runtime_error(errorMessage);
                                }
                            }
                        }
                    }
                }
            }
            catch(YAML::BadFile&e)
            {
                std::cout << "YAML Exception caught: " << e.what() << std::endl;
            }
        }

};

#endif // SYSTEMPARAMETERS_HPP
