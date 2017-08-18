#ifndef OBSTACLE_HPP
#define OBSTACLE_HPP

#include <iostream>
#include <armadillo>

#include <yaml-cpp/yaml.h>

/**
 * @brief The Obstacle class
 */
template <uword obsDim>
class Obstacle
{
    public:

        /**
         * @brief Obstacle
         */
        Obstacle(): dim(-1)
        {
            pos.zeros();
            radius  = 1.0;
        }

        /**
         * @brief Obstacle
         * @param other
         */
        Obstacle(const Obstacle<obsDim>&other)
        {
            copy(other);
        }

        /**
         * @brief getDim
         * @return
         */
        int Dim() const
        {
            return dim;
        }

        /**
         * @brief copy
         * @param other
         */
        void copy(const Obstacle&other)
        {
            dim     = other.dim;
            pos     = other.pos;
            radius  = other.radius;
        }

        /**
     * @brief Get the x-coordinate
     * @return
     */
        double x() const
        {
            return pos(0);
        }

        /**
     * @brief Get the y-coordinate
     * @return
     */
        double y() const
        {
            return pos(1);
        }

        /**
         * @brief z
         * @return
         */
        double z() const
        {
            return pos(2);
        }

        /**
         * @brief getRadius
         * @return
         */
        double getRadius() const
        {
            return radius;
        }

    public:
        int  dim;
        arma::vec::fixed<obsDim> pos;
        double radius;
};


/**
 * @brief loadMapYAML
 * @param filename
 * @param obstacles
 * @param bottomLeft
 * @param topRight
 */
template<uword obsDim>
void loadMapYAML(const std::string&filename,
                 std::vector< Obstacle<obsDim> >&obstacles,
                 vec::fixed<obsDim>&bottomLeft,
                 vec::fixed<obsDim>&topRight)
{
    obstacles.clear();

    Obstacle<obsDim> obstacle;

    try
    {

        YAML::Node parentNode = YAML::LoadFile(filename);

        for (auto it_parent = parentNode.begin(); it_parent != parentNode.end(); it_parent++)
        {
            YAML::Node key      = it_parent->first;
            YAML::Node child    = it_parent->second;

            const std::string str_key   = key.as<std::string>();

            //std::cout << "Parent is " << str_key << "\n"; // Parent

            if (parentNode.IsMap())
            {

                if(child.IsMap())
                {

                    // Now that you've got a map, you can iterate through it:
                    for (auto it = child.begin(); it != child.end(); ++it)
                    {


                        const YAML::Node childKey     = it->first;
                        const YAML::Node childValue   = it->second;

                        if(str_key == "bottomLeft")
                        {
                            for(unsigned int i=0; i<obsDim; i++)
                            {
                                bottomLeft(i)   = childValue[i].as<double>();
                            }

                        }

                        else if(str_key == "topRight")
                        {
                            for(unsigned int i=0; i<obsDim; i++)
                            {
                                topRight(i)   = childValue[i].as<double>();
                            }

                        }

                        else
                        {

                            if(childKey.as<std::string>() == "position")
                            {

                                for(unsigned int i=0; i<obsDim; i++)
                                {
                                    obstacle.pos(i) = childValue[i].as<double>();
                                }

                            }

                            else if (childKey.as<std::string>() == "radius")
                            {
                                obstacle.radius = childValue.as<double>();

                            }

                            else if (childKey.as<std::string>() == "dim")
                            {
                                obstacle.dim = childValue.as<double>();

                            }

                        }

                    }

                    if(str_key.find("Obstacle") != std::string::npos)
                    {
                        obstacles.push_back(obstacle);
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

#endif // OBSTACLE_HPP
