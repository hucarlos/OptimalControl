/****************************************************************************
**
** Copyright (C) 2014 Klaralvdalens Datakonsult AB (KDAB).
** Contact: https://www.qt.io/licensing/
**
** This file is part of the Qt3D module of the Qt Toolkit.
**
** $QT_BEGIN_LICENSE:BSD$
** Commercial License Usage
** Licensees holding valid commercial Qt licenses may use this file in
** accordance with the commercial license agreement provided with the
** Software or, alternatively, in accordance with the terms contained in
** a written agreement between you and The Qt Company. For licensing terms
** and conditions see https://www.qt.io/terms-conditions. For further
** information use the contact form at https://www.qt.io/contact-us.
**
** BSD License Usage
** Alternatively, you may use this file under the terms of the BSD license
** as follows:
**
** "Redistribution and use in source and binary forms, with or without
** modification, are permitted provided that the following conditions are
** met:
**   * Redistributions of source code must retain the above copyright
**     notice, this list of conditions and the following disclaimer.
**   * Redistributions in binary form must reproduce the above copyright
**     notice, this list of conditions and the following disclaimer in
**     the documentation and/or other materials provided with the
**     distribution.
**   * Neither the name of The Qt Company Ltd nor the names of its
**     contributors may be used to endorse or promote products derived
**     from this software without specific prior written permission.
**
**
** THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
** "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
** LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
** A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
** OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
** SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
** LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
** DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
** THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
** (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
** OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE."
**
** $QT_END_LICENSE$
**
****************************************************************************/

#include "scenemodifier.h"

#include <QtCore/QDebug>
#include <iostream>
#include <Utils/Utils.hpp>

#include <QSceneLoader>
#include <QtMath>

SceneModifier::SceneModifier(Qt3DCore::QEntity *rootEntity)
    : m_rootEntity(rootEntity)
{

    int redRGB  = 0xCC0000;
    robotColor  = QColor(QRgb(0xCCCC00));
    isReg       = true;

    // Cylinder shape data
    Qt3DExtras::QCylinderMesh *cylinder1 = new Qt3DExtras::QCylinderMesh();
    cylinder1->setRadius(0.4);
    cylinder1->setLength(6);
    cylinder1->setRings(100);
    cylinder1->setSlices(20);

    // CylinderMesh Transform
    Qt3DCore::QTransform *cylinderTransform1 = new Qt3DCore::QTransform();
    cylinderTransform1->setScale(1.0f);
    cylinderTransform1->setRotation(QQuaternion::fromAxisAndAngle(QVector3D(1.0f, 0.0f, 0.0f), 90.0f));
    cylinderTransform1->setTranslation(QVector3D(-1.25, 1.25, 0.0));

    Qt3DExtras::QPhongMaterial *cylinderMaterial1 = new Qt3DExtras::QPhongMaterial();
    cylinderMaterial1->setDiffuse(QColor(QRgb(redRGB)));

    // Cylinder
    m_cylinderEntity1 = new Qt3DCore::QEntity(m_rootEntity);
    m_cylinderEntity1->addComponent(cylinder1);
    m_cylinderEntity1->addComponent(cylinderMaterial1);
    m_cylinderEntity1->addComponent(cylinderTransform1);

    // ==========================================================================================
    Qt3DExtras::QCylinderMesh *cylinder2 = new Qt3DExtras::QCylinderMesh();
    cylinder2->setRadius(0.4);
    cylinder2->setLength(6);
    cylinder2->setRings(100);
    cylinder2->setSlices(20);

    // CylinderMesh Transform
    Qt3DCore::QTransform *cylinderTransform2 = new Qt3DCore::QTransform();
    cylinderTransform2->setScale(1.0f);
    cylinderTransform2->setRotation(QQuaternion::fromAxisAndAngle(QVector3D(1.0f, 0.0f, 0.0f), 90.0f));
    cylinderTransform2->setTranslation(QVector3D(1.25, 1.25, 0.0));

    Qt3DExtras::QPhongMaterial *cylinderMaterial2 = new Qt3DExtras::QPhongMaterial();
    cylinderMaterial2->setDiffuse(QColor(QRgb(redRGB)));

    // Cylinder
    m_cylinderEntity2 = new Qt3DCore::QEntity(m_rootEntity);
    m_cylinderEntity2->addComponent(cylinder2);
    m_cylinderEntity2->addComponent(cylinderMaterial2);
    m_cylinderEntity2->addComponent(cylinderTransform2);

    // ==========================================================================================
    Qt3DExtras::QCylinderMesh *cylinder3 = new Qt3DExtras::QCylinderMesh();
    cylinder3->setRadius(0.4);
    cylinder3->setLength(6);
    cylinder3->setRings(100);
    cylinder3->setSlices(20);

    // CylinderMesh Transform
    Qt3DCore::QTransform *cylinderTransform3 = new Qt3DCore::QTransform();
    cylinderTransform3->setScale(1.0f);
    cylinderTransform3->setRotation(QQuaternion::fromAxisAndAngle(QVector3D(1.0f, 0.0f, 0.0f), 90.0f));
    cylinderTransform3->setTranslation(QVector3D(1.25, -1.25, 0.0));

    Qt3DExtras::QPhongMaterial *cylinderMaterial3 = new Qt3DExtras::QPhongMaterial();
    cylinderMaterial3->setDiffuse(QColor(QRgb(redRGB)));

    // Cylinder
    m_cylinderEntity3 = new Qt3DCore::QEntity(m_rootEntity);
    m_cylinderEntity3->addComponent(cylinder3);
    m_cylinderEntity3->addComponent(cylinderMaterial3);
    m_cylinderEntity3->addComponent(cylinderTransform3);


    // ==========================================================================================
    Qt3DExtras::QCylinderMesh *cylinder4 = new Qt3DExtras::QCylinderMesh();
    cylinder4->setRadius(0.4);
    cylinder4->setLength(6);
    cylinder4->setRings(100);
    cylinder4->setSlices(20);

    // CylinderMesh Transform
    Qt3DCore::QTransform *cylinderTransform4 = new Qt3DCore::QTransform();
    cylinderTransform4->setScale(1.0f);
    cylinderTransform4->setRotation(QQuaternion::fromAxisAndAngle(QVector3D(1.0f, 0.0f, 0.0f), 90.0f));
    cylinderTransform4->setTranslation(QVector3D(-1.25, -1.25, 0.0));

    Qt3DExtras::QPhongMaterial *cylinderMaterial4 = new Qt3DExtras::QPhongMaterial();
    cylinderMaterial4->setDiffuse(QColor(QRgb(redRGB)));

    // Cylinder
    m_cylinderEntity4 = new Qt3DCore::QEntity(m_rootEntity);
    m_cylinderEntity4->addComponent(cylinder4);
    m_cylinderEntity4->addComponent(cylinderMaterial4);
    m_cylinderEntity4->addComponent(cylinderTransform4);


    // ==========================================================================================
    Qt3DExtras::QCylinderMesh *cylinder5 = new Qt3DExtras::QCylinderMesh();
    cylinder5->setRadius(0.4);
    cylinder5->setLength(6);
    cylinder5->setRings(100);
    cylinder5->setSlices(20);

    // CylinderMesh Transform
    Qt3DCore::QTransform *cylinderTransform5 = new Qt3DCore::QTransform();
    cylinderTransform5->setScale(1.0f);
    cylinderTransform5->setRotation(QQuaternion::fromAxisAndAngle(QVector3D(1.0f, 0.0f, 0.0f), 0.0f));
    cylinderTransform5->setTranslation(QVector3D(-1.25, 0.0, 1.25));

    Qt3DExtras::QPhongMaterial *cylinderMaterial5 = new Qt3DExtras::QPhongMaterial();
    cylinderMaterial5->setDiffuse(QColor(QRgb(redRGB)));

    // Cylinder
    m_cylinderEntity5 = new Qt3DCore::QEntity(m_rootEntity);
    m_cylinderEntity5->addComponent(cylinder5);
    m_cylinderEntity5->addComponent(cylinderMaterial5);
    m_cylinderEntity5->addComponent(cylinderTransform5);

    // ==========================================================================================
    Qt3DExtras::QCylinderMesh *cylinder6 = new Qt3DExtras::QCylinderMesh();
    cylinder6->setRadius(0.4);
    cylinder6->setLength(6);
    cylinder6->setRings(100);
    cylinder6->setSlices(20);

    // CylinderMesh Transform
    Qt3DCore::QTransform *cylinderTransform6 = new Qt3DCore::QTransform();
    cylinderTransform6->setScale(1.0f);
    cylinderTransform6->setRotation(QQuaternion::fromAxisAndAngle(QVector3D(1.0f, 0.0f, 0.0f), 0.0f));
    cylinderTransform6->setTranslation(QVector3D(1.25, 0.0, 1.25));

    Qt3DExtras::QPhongMaterial *cylinderMaterial6 = new Qt3DExtras::QPhongMaterial();
    cylinderMaterial6->setDiffuse(QColor(QRgb(redRGB)));

    // Cylinder
    m_cylinderEntity6 = new Qt3DCore::QEntity(m_rootEntity);
    m_cylinderEntity6->addComponent(cylinder6);
    m_cylinderEntity6->addComponent(cylinderMaterial6);
    m_cylinderEntity6->addComponent(cylinderTransform6);

    // ==========================================================================================
    Qt3DExtras::QCylinderMesh *cylinder7 = new Qt3DExtras::QCylinderMesh();
    cylinder7->setRadius(0.4);
    cylinder7->setLength(6);
    cylinder7->setRings(100);
    cylinder7->setSlices(20);

    // CylinderMesh Transform
    Qt3DCore::QTransform *cylinderTransform7 = new Qt3DCore::QTransform();
    cylinderTransform7->setScale(1.0f);
    cylinderTransform7->setRotation(QQuaternion::fromAxisAndAngle(QVector3D(1.0f, 0.0f, 0.0f), 0.0f));
    cylinderTransform7->setTranslation(QVector3D(1.25, 0.0, -1.25));

    Qt3DExtras::QPhongMaterial *cylinderMaterial7 = new Qt3DExtras::QPhongMaterial();
    cylinderMaterial7->setDiffuse(QColor(QRgb(redRGB)));

    // Cylinder
    m_cylinderEntity7 = new Qt3DCore::QEntity(m_rootEntity);
    m_cylinderEntity7->addComponent(cylinder7);
    m_cylinderEntity7->addComponent(cylinderMaterial7);
    m_cylinderEntity7->addComponent(cylinderTransform7);

    // ==========================================================================================
    Qt3DExtras::QCylinderMesh *cylinder8 = new Qt3DExtras::QCylinderMesh();
    cylinder8->setRadius(0.4);
    cylinder8->setLength(6);
    cylinder8->setRings(100);
    cylinder8->setSlices(20);

    // CylinderMesh Transform
    Qt3DCore::QTransform *cylinderTransform8 = new Qt3DCore::QTransform();
    cylinderTransform8->setScale(1.0f);
    cylinderTransform8->setRotation(QQuaternion::fromAxisAndAngle(QVector3D(1.0f, 0.0f, 0.0f), 0.0f));
    cylinderTransform8->setTranslation(QVector3D(-1.25, 0.0, -1.25));

    Qt3DExtras::QPhongMaterial *cylinderMaterial8 = new Qt3DExtras::QPhongMaterial();
    cylinderMaterial8->setDiffuse(QColor(QRgb(0xCC0000)));

    // Cylinder
    m_cylinderEntity8 = new Qt3DCore::QEntity(m_rootEntity);
    m_cylinderEntity8->addComponent(cylinder8);
    m_cylinderEntity8->addComponent(cylinderMaterial8);
    m_cylinderEntity8->addComponent(cylinderTransform8);

    // ==========================================================================================
    Qt3DExtras::QCylinderMesh *cylinder9 = new Qt3DExtras::QCylinderMesh();
    cylinder9->setRadius(0.4);
    cylinder9->setLength(6);
    cylinder9->setRings(100);
    cylinder9->setSlices(20);

    // CylinderMesh Transform
    Qt3DCore::QTransform *cylinderTransform9 = new Qt3DCore::QTransform();
    cylinderTransform9->setScale(1.0f);
    cylinderTransform9->setRotation(QQuaternion::fromAxisAndAngle(QVector3D(0.0f, 0.0f, 1.0f), 90.0f));
    cylinderTransform9->setTranslation(QVector3D(0.0, -1.25, 1.25));

    Qt3DExtras::QPhongMaterial *cylinderMaterial9 = new Qt3DExtras::QPhongMaterial();
    cylinderMaterial9->setDiffuse(QColor(QRgb(redRGB)));

    // Cylinder
    m_cylinderEntity9 = new Qt3DCore::QEntity(m_rootEntity);
    m_cylinderEntity9->addComponent(cylinder9);
    m_cylinderEntity9->addComponent(cylinderMaterial9);
    m_cylinderEntity9->addComponent(cylinderTransform9);

    // ==========================================================================================
    Qt3DExtras::QCylinderMesh *cylinder10 = new Qt3DExtras::QCylinderMesh();
    cylinder10->setRadius(0.4);
    cylinder10->setLength(6);
    cylinder10->setRings(100);
    cylinder10->setSlices(20);

    // CylinderMesh Transform
    Qt3DCore::QTransform *cylinderTransform10 = new Qt3DCore::QTransform();
    cylinderTransform10->setScale(1.0f);
    cylinderTransform10->setRotation(QQuaternion::fromAxisAndAngle(QVector3D(0.0f, 0.0f, 1.0f), 90.0f));
    cylinderTransform10->setTranslation(QVector3D(0.0, 1.25, 1.25));

    Qt3DExtras::QPhongMaterial *cylinderMaterial10 = new Qt3DExtras::QPhongMaterial();
    cylinderMaterial10->setDiffuse(QColor(QRgb(redRGB)));

    // Cylinder
    m_cylinderEntity10 = new Qt3DCore::QEntity(m_rootEntity);
    m_cylinderEntity10->addComponent(cylinder10);
    m_cylinderEntity10->addComponent(cylinderMaterial10);
    m_cylinderEntity10->addComponent(cylinderTransform10);

    // ==========================================================================================
    Qt3DExtras::QCylinderMesh *cylinder11 = new Qt3DExtras::QCylinderMesh();
    cylinder11->setRadius(0.4);
    cylinder11->setLength(6);
    cylinder11->setRings(100);
    cylinder11->setSlices(20);

    // CylinderMesh Transform
    Qt3DCore::QTransform *cylinderTransform11 = new Qt3DCore::QTransform();
    cylinderTransform11->setScale(1.0f);
    cylinderTransform11->setRotation(QQuaternion::fromAxisAndAngle(QVector3D(0.0f, 0.0f, 1.0f), 90.0f));
    cylinderTransform11->setTranslation(QVector3D(0.0, 1.25, -1.25));

    Qt3DExtras::QPhongMaterial *cylinderMaterial11 = new Qt3DExtras::QPhongMaterial();
    cylinderMaterial11->setDiffuse(QColor(QRgb(redRGB)));

    // Cylinder
    m_cylinderEntity11 = new Qt3DCore::QEntity(m_rootEntity);
    m_cylinderEntity11->addComponent(cylinder11);
    m_cylinderEntity11->addComponent(cylinderMaterial11);
    m_cylinderEntity11->addComponent(cylinderTransform11);

    // ==========================================================================================
    Qt3DExtras::QCylinderMesh *cylinder12 = new Qt3DExtras::QCylinderMesh();
    cylinder12->setRadius(0.4);
    cylinder12->setLength(6);
    cylinder12->setRings(100);
    cylinder12->setSlices(20);

    // CylinderMesh Transform
    Qt3DCore::QTransform *cylinderTransform12 = new Qt3DCore::QTransform();
    cylinderTransform12->setScale(1.0f);
    cylinderTransform12->setRotation(QQuaternion::fromAxisAndAngle(QVector3D(0.0f, 0.0f, 1.0f), 90.0f));
    cylinderTransform12->setTranslation(QVector3D(0.0, -1.25, -1.25));

    Qt3DExtras::QPhongMaterial *cylinderMaterial12 = new Qt3DExtras::QPhongMaterial();
    cylinderMaterial12->setDiffuse(QColor(QRgb(redRGB)));

    // Cylinder
    m_cylinderEntity12 = new Qt3DCore::QEntity(m_rootEntity);
    m_cylinderEntity12->addComponent(cylinder12);
    m_cylinderEntity12->addComponent(cylinderMaterial12);
    m_cylinderEntity12->addComponent(cylinderTransform12);


    //Robot
    robot  = new Qt3DCore::QEntity(m_rootEntity);
    Qt3DExtras::QSphereMesh *sphereMesh = new Qt3DExtras::QSphereMesh;
    sphereMesh->setRadius(0.3);

    robotMaterial = new Qt3DExtras::QPhongMaterial();
    robotMaterial->setDiffuse(robotColor);


    robotTransform  = new Qt3DCore::QTransform;
    controller      = new RobotController(robotTransform);
    controller->setTarget(robotTransform);
    controller->setRadius(1.0f);
    controller->updateMatrix();

//    robot->addComponent(sphereMesh);
//    robot->addComponent(robotTransform);
//    robot->addComponent(robotMaterial);

    //=================================================== SCENE==============================



    Qt3DCore::QEntity *sceneEntity = new Qt3DCore::QEntity(m_rootEntity);
    Qt3DRender::QSceneLoader *scene = new Qt3DRender::QSceneLoader();
    scene->setObjectName(QStringLiteral("scene"));
    QUrl url = QUrl::fromLocalFile("Models/quadrotorsmall/quadrotor1.obj");
    scene->setSource(url);

    sceneTransform = new Qt3DCore::QTransform();

    sceneEntity->addComponent(sceneTransform);

    sceneEntity->addComponent(scene);


    //Load data
    data1.load("SELQRarma2.txt",    arma::raw_ascii);
    data2.load("iQRSELQRarma2.txt", arma::raw_ascii);

    k=0;

    timer = new QTimer(this);
    connect(timer, SIGNAL(timeout()), this, SLOT(update()));
    timer->start(40);


}

void SceneModifier::update()
{
    if(k < data1.n_rows)
    {
        if(isReg)
        {
            controller->setX(data2(k,0));
            controller->setY(data2(k,1));
            controller->setZ(data2(k,2));

            sceneTransform->setScale(0.001);
            sceneTransform->setTranslation(QVector3D(data2(k,0), data2(k,1), data2(k,2)));
            sceneTransform->setRotation(QQuaternion::fromAxisAndAngle(QVector3D(1.0f, 0.0f, 0.0f), qRadiansToDegrees(data2(k,6))));
            sceneTransform->setRotation(QQuaternion::fromAxisAndAngle(QVector3D(0.0f, 1.0f, 0.0f), qRadiansToDegrees(data2(k,7))));
            sceneTransform->setRotation(QQuaternion::fromAxisAndAngle(QVector3D(0.0f, 0.0f, 1.0f), qRadiansToDegrees(data2(k,8))));
        }

        else
        {
            controller->setX(data1(k,0));
            controller->setY(data1(k,1));
            controller->setZ(data1(k,2));

            sceneTransform->setScale(0.001);
            sceneTransform->setTranslation(QVector3D(data1(k,0), data1(k,1), data1(k,2)));
            sceneTransform->setRotation(QQuaternion::fromAxisAndAngle(QVector3D(1.0f, 0.0f, 0.0f), qRadiansToDegrees(data1(k,6))));
            sceneTransform->setRotation(QQuaternion::fromAxisAndAngle(QVector3D(0.0f, 1.0f, 0.0f), qRadiansToDegrees(data1(k,7))));
            sceneTransform->setRotation(QQuaternion::fromAxisAndAngle(QVector3D(0.0f, 0.0f, 1.0f), qRadiansToDegrees(data1(k,8))));
        }
        controller->updateMatrix();
        k++;
    }
    else
    {
        k=0;
    }

}

SceneModifier::~SceneModifier()
{
}

/**
 * @brief SceneModifier::enableRegression
 * @param enabled
 */
void SceneModifier::enableRegression(bool enabled)
{
    isReg       = enabled;
    k           = 0;
    if(enabled)
        robotColor = QColor(QRgb(0xCCCC00));
    else
        robotColor = QColor(QRgb(0x004C99));

    robotMaterial->setDiffuse(robotColor);
}
