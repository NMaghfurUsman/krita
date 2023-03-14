/*
 * SPDX-FileCopyrightText: 2008 Cyrille Berger <cberger@cberger.net>
 * SPDX-FileCopyrightText: 2010 Geoffry Song <goffrie@gmail.com>
 * SPDX-FileCopyrightText: 2021 Nabil Maghfur Usman <nmaghfurusman@gmail.com>
 *
 *  SPDX-License-Identifier: LGPL-2.0-or-later
 */

#include "ThreePointAssistant.h"
#include "kis_debug.h"
#include "kis_painting_assistant.h"
#include <klocalizedstring.h>

#include <QPainter>
#include <QPainterPath>
#include <QLinearGradient>
#include <QTransform>

#include <kis_canvas2.h>
#include <kis_coordinates_converter.h>
#include <kis_algebra_2d.h>
#include <kis_dom_utils.h>
#include <math.h>
#include <QtCore/qmath.h>
#include <kis_assert.h>
#include <QtMath>
#include <qglobal.h>
#include <QVector3D>
#include <QQuaternion>

class ThreePoint
{
    bool isValid();
    QPointF orthocenter();
    qreal radius();
    QTransform localTransform;
};

inline static QLineF altitude(QPointF a, QPointF b, QPointF c)
{
    return QLineF(a, c).normalVector().translated(b-a);
}

inline static QPointF orthocenter(QPointF a, QPointF b, QPointF c)
{
    QLineF altitudeA = altitude(a,b,c);
    QLineF altitudeB = altitude(b,a,c);
    QPointF intersection;
    altitudeA.intersect(altitudeB, &intersection);
    return intersection;
}

ThreePointAssistant::ThreePointAssistant()
    : NPointPerspective("three point", i18n("Three point assistant"))
{
}

ThreePointAssistant::ThreePointAssistant(const ThreePointAssistant &rhs, QMap<KisPaintingAssistantHandleSP, KisPaintingAssistantHandleSP> &handleMap)
    : NPointPerspective(rhs, handleMap)
    , m_canvas(rhs.m_canvas)
    , m_snapLine(rhs.m_snapLine)
    , m_gridDensity(rhs.m_gridDensity)
    , m_useVertical(rhs.m_useVertical)
    , m_lastUsedPoint(rhs.m_lastUsedPoint)
{
}

KisPaintingAssistantSP ThreePointAssistant::clone(QMap<KisPaintingAssistantHandleSP, KisPaintingAssistantHandleSP> &handleMap) const
{
    return KisPaintingAssistantSP(new ThreePointAssistant(*this, handleMap));
}

QPointF ThreePointAssistant::project(const QPointF& point, const QPointF& strokeBegin, const bool snapToAny)
{
    Q_ASSERT(isAssistantComplete());

    return point;
}

void ThreePointAssistant::endStroke()
{
    m_lastUsedPoint = -1;
    KisPaintingAssistant::endStroke();
}

QPointF ThreePointAssistant::adjustPosition(const QPointF& pt, const QPointF& strokeBegin, const bool snapToAny)
{
    return project(pt, strokeBegin, snapToAny);
}

void ThreePointAssistant::adjustLine(QPointF &point, QPointF &strokeBegin)
{
    QPointF p = project(point, strokeBegin, true);
    point = p;
}

void ThreePointAssistant::drawAssistant(QPainter& gc, const QRectF& updateRect, const KisCoordinatesConverter* converter, bool cached, KisCanvas2* canvas, bool assistantVisible, bool previewVisible)
{
    Q_UNUSED(updateRect);
    Q_UNUSED(cached);
    gc.save();
    gc.resetTransform();
    QPointF mousePos(0,0);

    const QTransform initialTransform = converter->documentToWidgetTransform();
    bool isEditing = false;
    bool showLocal = isLocal() && handles().size() == 5;

    if (canvas){
        //simplest, cheapest way to get the mouse-position//
        mousePos = canvas->canvasWidget()->mapFromGlobal(QCursor::pos());
        isEditing = canvas->paintingAssistantsDecoration()->isEditingAssistants();
        m_canvas = canvas;
    }
    else {
        mousePos = QCursor::pos();//this'll give an offset//
        dbgFile<<"canvas does not exist in ruler, you may have passed arguments incorrectly:"<<canvas;
    }

    if (m_followBrushPosition && m_adjustedPositionValid) {
        mousePos = initialTransform.map(m_adjustedBrushPosition);
    }

    if (isEditing) {
        Q_FOREACH (const QPointF* handle, handles()) {
            QPointF h = initialTransform.map(*handle);
            QRectF ellipse = QRectF(QPointF(h.x() -15, h.y() -15), QSizeF(30, 30));

            QPainterPath pathCenter;
            pathCenter.addEllipse(ellipse);
            drawPath(gc, pathCenter, isSnappingActive());
        }
    }

    if (handles().size() >= 2) {
        const QPointF p1 = *handles()[0];
        const QPointF p2 = *handles()[1];
        const QRect viewport= gc.viewport();

        const QPolygonF localPoly = (isLocal() && handles().size() == 5) ? initialTransform.map(QPolygonF(getLocalRect())) : QPolygonF();
        const QPolygonF viewportAndLocalPoly = !localPoly.isEmpty() ? QPolygonF(QRectF(viewport)).intersected(localPoly) : QRectF(viewport);


        QPainterPath path;
        QPainterPath previewPath; // part of the preview, instead of the assistant itself

        // draw the horizon
        if (assistantVisible == true || isEditing == true) {
            QLineF horizonLine = initialTransform.map(QLineF(p1,p2));
            KisAlgebra2D::cropLineToConvexPolygon(horizonLine, viewportAndLocalPoly, true, true);
            path.moveTo(horizonLine.p1());
            path.lineTo(horizonLine.p2());
        }

        // draw the VP-->mousePos lines
        if (isEditing == false && previewVisible == true && isSnappingActive() == true) {
            // draw the line vp <-> mouse even outside of the local rectangle
            // but only if the mouse pos is inside the rectangle
            QLineF snapMouse1 = QLineF(initialTransform.map(p1), mousePos);
            QLineF snapMouse2 = QLineF(initialTransform.map(p2), mousePos);
            KisAlgebra2D::cropLineToConvexPolygon(snapMouse1, viewportAndLocalPoly, false, true);
            KisAlgebra2D::cropLineToConvexPolygon(snapMouse2, viewportAndLocalPoly, false, true);
            previewPath.moveTo(snapMouse1.p1());
            previewPath.lineTo(snapMouse1.p2());
            previewPath.moveTo(snapMouse2.p1());
            previewPath.lineTo(snapMouse2.p2());
        }

        // draw the local rectangle
        if (showLocal && assistantVisible) {
            QPointF p1 = *handles()[(int)LocalFirstHandle];
            QPointF p3 = *handles()[(int)LocalSecondHandle];
            QPointF p2 = QPointF(p1.x(), p3.y());
            QPointF p4 = QPointF(p3.x(), p1.y());

            path.moveTo(initialTransform.map(p1));

            path.lineTo(initialTransform.map(p2));
            path.lineTo(initialTransform.map(p3));
            path.lineTo(initialTransform.map(p4));
            path.lineTo(initialTransform.map(p1));
        }
        // draw temp vanishing point lines
        if (!isValid()) {
            QPainterPath tempPath;
            int tempDensity = m_gridDensity * 20; // the vanishing point density seems visibly more dense, hence let's make it less dense
            for (int i = 0; i < handles().size() && i < 3; i++) { // only draw on first 3 handles
                const QPointF p = initialTransform.map(*handles()[i]);
                for (int currentAngle=0; currentAngle <= 180; currentAngle = currentAngle + tempDensity) {

                    // determine the correct angle based on the iteration
                    float xPos = cos(currentAngle * M_PI / 180);
                    float yPos = sin(currentAngle * M_PI / 180);
                    float length = 100;
                    QPointF unit = QPointF(length*xPos, length*yPos);

                    // find point
                    QLineF vpLine = QLineF(p - unit, p + unit);
                    if (KisAlgebra2D::intersectLineRect(vpLine, viewport, false)) {
                        // make a line from VP center to edge of canvas with that angle
                        tempPath.moveTo(vpLine.p1());
                        tempPath.lineTo(vpLine.p2());
                    }
                }
                drawPath(gc, tempPath, true, true);
            }
        }


        if (handles().size() >= 3) {
            const QPointF p3 = *handles()[2];

            if (assistantVisible && isEditing) {
                drawX(gc, initialTransform.map(p1));
                drawX(gc, initialTransform.map(p2));
                drawX(gc, initialTransform.map(p3));
            }

            if (isValid() && assistantVisible) {

                const QTransform t = localTransform(p1,p2,p3);
                const QPointF vp_a = t.map(p1);
                const QPointF vp_b = t.map(p2);

                const QTransform inv = t.inverted()*initialTransform;

                const QPointF ortho = orthocenter(vp_a, vp_b, t.map(p3));

                const qreal mid = QLineF(vp_a, vp_b).pointAt(0.5).x();
                const qreal radius = QLineF(vp_a, vp_b).length() / 2.0;
                const qreal sp_distance = sqrt(radius*radius - mid*mid);
                const qreal theta = vp_a.y() - ortho.y();
                const qreal cov_size = sqrt(sp_distance*sp_distance - theta*theta);

                if (assistantVisible && isEditing) {
                    drawX(gc, inv.map(ortho)); // center of vision
                    const qreal actual_size = inv.map(QLineF(ortho, QPointF(cov_size, ortho.y()))).length();
                    const QPointF cov_corner = QPointF(actual_size,actual_size);
                    const QRectF cov_circle = QRectF(inv.map(ortho) - cov_corner, inv.map(ortho) + cov_corner);
                    path.moveTo(inv.map(ortho));
                    path.addEllipse(cov_circle); // 90 degree cone of vision
                }

                drawPath(gc, path, isSnappingActive(), true);
                path = QPainterPath();

                // orthocenter == image center
                // orthocenter as origin makes calculations easier to reason with
                QTransform t_ortho = QTransform();
                t_ortho.translate(-ortho.x(), -ortho.y());

                const QPointF principal_pt = QPointF(0,vp_a.y());
                const qreal principal_pt_dst = ortho.y() - principal_pt.y();

                // distance to image plane
                const qreal dst = principal_pt_dst / qTan(qAsin(principal_pt_dst / sp_distance));

                // orthonormal vectors projecting from center of projection to VPs on image plane
                const QVector3D p1_vec = QVector3D(QVector2D(t_ortho.map(vp_a)),dst);
                const QVector3D p3_vec = QVector3D(QVector2D(t_ortho.map(QPointF(0,0))),dst);

                // use the vanishing point vectors to recover the camera's orientation
                const qreal pitch_angle = qRadiansToDegrees(qAsin(qAbs(p3_vec.z()/p3_vec.length())));
                const qreal yaw_angle = qRadiansToDegrees(qAsin(qAbs(p1_vec.x()/p1_vec.length())));
                const QQuaternion pitch = QQuaternion::fromAxisAndAngle(QVector3D(1,0,0),
                                                                        // hack alert: not sure why it works like this
                                                                        pitch_angle * (principal_pt.y() < 0 ? 1 : -1));
                const QQuaternion yaw = QQuaternion::fromAxisAndAngle(QVector3D(0,1,0), -yaw_angle);
                const QQuaternion orientation = pitch * yaw;

                // draw perspective grid as 3d cube room centered around viewer
                const qreal unit_size = dst; // size of grid square
                const int plane_dst = 5; // distance from center of room to walls in unit size
                const QVector3D offset = pitch.rotatedVector(QVector3D(0,0,unit_size)); // offset from center of rotation
                const int wall_size = plane_dst; // wall size as units from center (so double this for actual size)

                for (int i = -wall_size; i <= wall_size; i++) {

                    for (int corner : {1,-1}) { // repeat for near and far corner of cube
                        // 6 edges = 1 half-cube corner
                        const std::array<QVector3D,2> edges[] = {
                            {QVector3D( corner*plane_dst , wall_size        , i                ),
                             QVector3D( corner*plane_dst ,-wall_size        , i                )},
                            {QVector3D( corner*plane_dst , i                , wall_size        ),
                             QVector3D( corner*plane_dst , i                ,-wall_size        )},
                            {QVector3D( wall_size        , i                , corner*plane_dst ),
                             QVector3D(-wall_size        , i                , corner*plane_dst )},
                            {QVector3D( i                , wall_size        , corner*plane_dst ),
                             QVector3D( i                ,-wall_size        , corner*plane_dst )},
                            {QVector3D( wall_size        , corner*plane_dst , i                ),
                             QVector3D(-wall_size        , corner*plane_dst , i                )},
                            {QVector3D( i                , corner*plane_dst , wall_size        ),
                             QVector3D( i                , corner*plane_dst ,-wall_size        )}
                        };

                        for (std::array<QVector3D,2> edge : edges) {
                            const QVector3D pt_far  = orientation.rotatedVector(edge[0] * unit_size) + offset;
                            const QVector3D pt_near = orientation.rotatedVector(edge[1] * unit_size) + offset;
                            if (pt_far.z() > 0 || pt_near.z() > 0) { // don't draw line if it's behind image plane
                                const QPointF projection_far  = QPointF(pt_far.x()/pt_far.z()*dst, pt_far.y()/pt_far.z()*dst);
                                const QPointF projection_near = QPointF(pt_near.x()/pt_near.z()*dst, pt_near.y()/pt_near.z()*dst);
                                const QPointF projection_far_mapped  = inv.map(t_ortho.inverted().map(projection_far));
                                const QPointF projection_near_mapped = inv.map(t_ortho.inverted().map(projection_near));
                                QLineF grid_line = QLineF(projection_near_mapped, projection_far_mapped);
                                KisAlgebra2D::cropLineToConvexPolygon(grid_line, viewportAndLocalPoly,
                                                                      // extend line in opposite direction if pt is behind image plane
                                                                      pt_far.z() < 0, pt_near.z() < 0);
                                if (!grid_line.isNull()) {
                                    // draw in opposite direction if one of the points is behind the image plane
                                    if (pt_far.z() < 0 || pt_near.z() < 0) {
                                        if (pt_near.z() < 0) {grid_line.setP1(projection_far_mapped);}
                                        if (pt_far.z() < 0) {grid_line.setP2(projection_near_mapped);}
                                        KisAlgebra2D::cropLineToConvexPolygon(grid_line, viewportAndLocalPoly, false, false);
                                    }

                                    path.moveTo(grid_line.p1());
                                    path.lineTo(grid_line.p2());
                                }
                            }
                        }
                    }
                }

                drawPath(gc, path, isSnappingActive(), true);
                path = QPainterPath();
            } else {
                drawError(gc, path);
            }
        } else {
            drawPreview(gc, path);
        }
    }

    gc.restore();
    //KisPaintingAssistant::drawAssistant(gc, updateRect, converter, cached, canvas, assistantVisible, previewVisible);
}


QTransform ThreePointAssistant::localTransform(QPointF vp_a, QPointF vp_b, QPointF pt_c)
{
    QTransform t = QTransform();
    t.rotate(QLineF(vp_a, vp_b).angle());
    t.translate(-pt_c.x(),-pt_c.y());
    return t;
}

bool ThreePointAssistant::isValid()
{
    if (handles().size() < 3) {
        return false;
    } else {
        const QTransform t = localTransform(*handles()[0],*handles()[1],*handles()[2]);
        const QPointF vp_a = t.map(*handles()[0]);
        const QPointF vp_b = t.map(*handles()[1]);

        bool vpMode = true;

        // is vertical vanishing point between the left and right vanishing points?
        bool isBetween = (vp_a.x() < 0 && vp_b.x() > 0) || (vp_a.x() > 0 && vp_b.x() < 0);

        // do the 3 handles form an acute triangle?
        QPointF vp_c = vpMode ? QPointF(0,0) : orthocenter(vp_a, vp_b, QPointF(0,0));
        QLineF a = QLineF(vp_c, vp_a);
        QLineF b = QLineF(vp_c, vp_b);
        bool isAcute = qCos(qDegreesToRadians(a.angleTo(b))) > 0;

        return isAcute && isBetween;
    }
}

void ThreePointAssistant::drawCache(QPainter& gc, const KisCoordinatesConverter *converter, bool assistantVisible)
{
    Q_UNUSED(gc);
    Q_UNUSED(converter);
    Q_UNUSED(assistantVisible);
    if (!m_canvas || !isAssistantComplete()) {
        return;
    }

    if (assistantVisible == false ||   m_canvas->paintingAssistantsDecoration()->isEditingAssistants()) {
        return;
    }
}

KisPaintingAssistantHandleSP ThreePointAssistant::firstLocalHandle() const
{
    if (handles().size() > LocalFirstHandle) {
        return handles().at(LocalFirstHandle);
    } else {
        return nullptr;
    }
}

KisPaintingAssistantHandleSP ThreePointAssistant::secondLocalHandle() const
{
    if (handles().size() > LocalSecondHandle) {
        return handles().at(LocalSecondHandle);
    } else {
        return nullptr;
    }
}

QPointF ThreePointAssistant::getDefaultEditorPosition() const
{
    if (!handles().empty()) {
        KIS_SAFE_ASSERT_RECOVER_RETURN_VALUE(false, *handles().at(0));
        return *handles().at(0);
    } else {
        KIS_SAFE_ASSERT_RECOVER_RETURN_VALUE(false, QPointF(0, 0));
        return QPointF(0, 0);
    }
}

void ThreePointAssistant::realignSideHandles(KisPaintingAssistantHandleSP dragged_handle) {
    Q_UNIMPLEMENTED();
}

void ThreePointAssistant::realignVanishingPoint(KisPaintingAssistantHandleSP dragged_handle, KoPointerEvent* event, QPointF* drag_start, QPointF* adjustment) {

    if (dragged_handle != handles()[VerticalHandle]) {
        KisPaintingAssistantHandleSP handleOpp =
            dragged_handle == handles()[FirstHandle] ? handles()[SecondHandle] : handles()[FirstHandle];
        const QPointF prevPoint = adjustment->isNull() ? *drag_start : *adjustment;

        const QTransform t = localTransform(prevPoint, *handleOpp, *handles()[VerticalHandle]);
        const QTransform inv = t.inverted();

        const QPointF vp_a = t.map(prevPoint);
        const QPointF vp_b = t.map(*handleOpp);
        const QLineF horizon = QLineF(vp_a, vp_b);
        const qreal mid = QLineF(vp_a, vp_b).pointAt(0.5).x();
        const qreal radius = QLineF(vp_a, vp_b).length() / 2.0;
        const qreal sp_distance = sqrt(radius*radius - mid*mid);
        const QPointF sp = QPointF(0, vp_a.y()+sp_distance);

        const bool preserve_distortion_snap = event->modifiers() == Qt::ControlModifier;
        QPointF snap_point;
        QPointF opp_snap_point;
        QLineF sp_to_opp_vp;

        if (preserve_distortion_snap) {
            const QLineF sp_to_vp = QLineF(sp, t.map(*dragged_handle));
            sp_to_opp_vp = sp_to_vp.normalVector();
            sp_to_vp.intersect(horizon,&snap_point);

            const bool no_intersection =
              // NB: opp_snap_point is initialized here
              sp_to_opp_vp.intersect(horizon, &opp_snap_point) == QLineF::NoIntersection;

            *dragged_handle = inv.map(snap_point);
            *handleOpp = inv.map(opp_snap_point);
        }
        *adjustment = *dragged_handle; // clear
    } else {
        const QPointF prevPoint = adjustment->isNull() ? *drag_start : *adjustment;

        const QTransform t = localTransform(*handles()[FirstHandle], *handles()[SecondHandle], prevPoint);
        const QTransform inv = t.inverted();

        const QPointF vp_a = t.map(*handles()[FirstHandle]);
        const QPointF vp_b = t.map(*handles()[SecondHandle]);
        const QLineF vertical = QLineF(QPointF(0,0), QPointF(0,vp_a.y()));
        const qreal mid = QLineF(vp_a, vp_b).pointAt(0.5).x();
        const qreal radius = QLineF(vp_a, vp_b).length() / 2.0;
        const qreal sp_distance = sqrt(radius*radius - mid*mid);
        const QPointF sp = QPointF(0, vp_a.y()+sp_distance);
        const QPointF ortho = orthocenter(vp_a, vp_b, t.map(prevPoint));
        const qreal theta = vp_a.y() - ortho.y();
        const qreal cov_size = sqrt(sp_distance*sp_distance - theta*theta);
        const bool preserve_distortion_snap = event->modifiers() == Qt::ControlModifier;

        QPointF snap_point;
        QLineF vsp_to_pp;
        QPointF pp;
        QPointF snap_point_a;
        QPointF snap_point_b;
        const QPointF vsp = QPointF(0,ortho.y()) + QPointF(cov_size,0);
        QLineF vsp_to_vvp = QLineF(vsp, t.map(*dragged_handle));

        if (preserve_distortion_snap) {
            vsp_to_vvp.intersect(vertical,&snap_point);
            snap_point = QPointF(0,snap_point.y());
            vsp_to_vvp = QLineF(vsp, snap_point);
            vsp_to_pp = vsp_to_vvp.normalVector();
            vsp_to_pp.intersect(vertical, &pp);
            vsp_to_pp.setP2(pp);
            const QLineF horizon = QLineF(pp, pp+QPointF(10,0));
            QLineF sp_to_vp_a = QLineF(sp, vp_a);
            QLineF sp_to_vp_b = QLineF(sp, vp_b);

            qreal new_sp_distance = vsp_to_pp.length();
            QPointF new_sp = QPointF(0, pp.y() + new_sp_distance);
            sp_to_vp_a.translate(new_sp - sp_to_vp_a.p1());
            sp_to_vp_b.translate(new_sp - sp_to_vp_b.p1());
            sp_to_vp_a.intersect(horizon, &snap_point_a);
            sp_to_vp_b.intersect(horizon, &snap_point_b);

            vsp_to_pp = QLineF(vsp,pp);
            vsp_to_vvp = vsp_to_pp.normalVector();
            vsp_to_vvp.intersect(vertical,&snap_point);

            const QLineF altitude_a = QLineF(snap_point_a, ortho);
            const QLineF altitude_b = QLineF(snap_point_b, ortho);
            QLineF vl_a = altitude_a.normalVector();
            vl_a.translate(snap_point_b - snap_point_a);
            QLineF vl_b = altitude_b.normalVector();
            vl_b.translate(snap_point_a - snap_point_b);
            vl_a.intersect(vl_b, &snap_point);

            *handles()[FirstHandle] = inv.map(snap_point_a);
            *handles()[SecondHandle] = inv.map(snap_point_b);
            *dragged_handle = inv.map(snap_point);
        }
        *adjustment = *dragged_handle; // clear
    }
}


void ThreePointAssistant::initSideHandles()
{
}

bool ThreePointAssistant::isAssistantComplete() const
{
    return handles().size() >= numHandles();
}

bool ThreePointAssistant::canBeLocal() const
{
    return true;
}
ThreePointAssistantFactory::ThreePointAssistantFactory()
{
}

ThreePointAssistantFactory::~ThreePointAssistantFactory()
{
}

QString ThreePointAssistantFactory::id() const
{
    return "three point";
}

QString ThreePointAssistantFactory::name() const
{
    return i18n("3 Point Perspective");
}

KisPaintingAssistant* ThreePointAssistantFactory::createPaintingAssistant() const
{
    return new ThreePointAssistant;
}
