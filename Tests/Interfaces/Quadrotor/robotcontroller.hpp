#ifndef ROBOTCONTROLLER_HPP
#define ROBOTCONTROLLER_HPP

#include <QObject>
#include <QMatrix4x4>

namespace Qt3DCore {
class QTransform;
}

class RobotController: public QObject
{
        Q_OBJECT
        Q_PROPERTY(Qt3DCore::QTransform* target READ target WRITE setTarget NOTIFY targetChanged)
        Q_PROPERTY(float radius READ radius WRITE setRadius NOTIFY radiusChanged)
        Q_PROPERTY(float angle READ angle WRITE setAngle NOTIFY angleChanged)

    public:
        RobotController();

        RobotController(QObject *parent = 0);

        void setTarget(Qt3DCore::QTransform *target);
        Qt3DCore::QTransform *target() const;

        void setRadius(float radius);
        float radius() const;

        void setAngle(float angle);
        float angle() const;

    signals:
        void targetChanged();
        void radiusChanged();
        void angleChanged();

    public:
        void updateMatrix();

        float x() const;
        void setX(float x);

        float y() const;

        float z() const;

        void setY(float y);

        void setZ(float z);

    private:
        Qt3DCore::QTransform *m_target;
        QMatrix4x4 m_matrix;
        float m_radius;
        float m_angle;

        // 3D positions
        float m_x;
        float m_y;
        float m_z;
};

#endif // ROBOTCONTROLLER_HPP
