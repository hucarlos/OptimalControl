#include "robotcontroller.hpp"

#include <Qt3DCore/qtransform.h>

QT_BEGIN_NAMESPACE

RobotController::RobotController(QObject *parent)
    : QObject(parent)
    , m_target(nullptr)
    , m_matrix()
    , m_radius(1.0f)
    , m_angle(0.0f)
{
}

void RobotController::setTarget(Qt3DCore::QTransform *target)
{
    if (m_target != target) {
        m_target = target;
        emit targetChanged();
    }
}

Qt3DCore::QTransform *RobotController::target() const
{
    return m_target;
}

void RobotController::setRadius(float radius)
{
    if (!qFuzzyCompare(radius, m_radius)) {
        m_x = radius;
        updateMatrix();
        emit radiusChanged();
    }
}

float RobotController::radius() const
{
    return m_radius;
}

void RobotController::setAngle(float angle)
{
    if (!qFuzzyCompare(angle, m_angle)) {
        m_angle = angle;
        updateMatrix();
        emit angleChanged();
    }
}

float RobotController::angle() const
{
    return m_angle;
}

void RobotController::updateMatrix()
{
    m_matrix.setToIdentity();
    m_matrix.rotate(m_angle, QVector3D(1.0f, 0.0f, 0.0f));
    m_matrix.translate(m_x, m_y, m_z);
    m_target->setMatrix(m_matrix);
}

float RobotController::x() const
{
    return m_x;
}

void RobotController::setX(float x)
{
    m_x = x;
}

float RobotController::y() const
{
    return m_y;
}

float RobotController::z() const
{
    return m_z;
}

void RobotController::setY(float y)
{
    m_y = y;
}

void RobotController::setZ(float z)
{
    m_z = z;
}

QT_END_NAMESPACE
