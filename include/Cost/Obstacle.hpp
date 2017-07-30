#ifndef OBSTACLE_HPP
#define OBSTACLE_HPP

#include <iostream>
#include <armadillo>

#include <yaml-cpp/yaml.h>

/**
 * @brief The Obstacle class
 */
template <int obsDim>
class Obstacle
{
    public:

        /**
         * @brief Obstacle
         */
        Obstacle(): dim(obsDim)
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
        size_t getDim() const
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
            return pos[0];
        }

        /**
     * @brief Get the y-coordinate
     * @return
     */
        double y() const
        {
            return pos[1];
        }

        double getRadius() const
        {
            return radius;
        }

    public:
        size_t dim;
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
template<int obsDim>
void loadMapYAML(const std::string&filename,
                 std::vector< Obstacle<obsDim> >&obstacles,
                 arma::vec2&bottomLeft,
                 arma::vec2&topRight)
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
                            bottomLeft(0)   = childValue[0].as<double>();
                            bottomLeft(1)   = childValue[1].as<double>();
                        }

                        else if(str_key == "topRight")
                        {
                            topRight(0)   = childValue[0].as<double>();
                            topRight(1)   = childValue[1].as<double>();
                        }

                        else
                        {

                            if(childKey.as<std::string>() == "position")
                            {

                                obstacle.pos(0) = childValue[0].as<double>();
                                obstacle.pos(1) = childValue[1].as<double>();

                            }

                            else if (childKey.as<std::string>() == "radius")
                            {
                                obstacle.radius = childValue.as<double>();

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
