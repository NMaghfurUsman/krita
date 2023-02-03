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
        const QRect viewport= gc.viewport();

        const QPolygonF localPoly = (isLocal() && handles().size() == 5) ? initialTransform.map(QPolygonF(getLocalRect())) : QPolygonF();
        const QPolygonF viewportAndLocalPoly = !localPoly.isEmpty() ? QPolygonF(QRectF(viewport)).intersected(localPoly) : QRectF(viewport);

        QPainterPath path;
        QPainterPath previewPath; // part of the preview, instead of the assistant itself
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

        // drawPreview(gc,previewPath);
        // drawPath(gc, path, isSnappingActive());

        if (handles().size() == 3) {
            // path = QPainterPath(); // clear
            const QPointF p3 = *handles()[2];

            if (assistantVisible && isEditing) {
                QLineF vanishingLineA = initialTransform.map(QLineF(p3,p1));
                QLineF vanishingLineB = initialTransform.map(QLineF(p3,p2));
                KisAlgebra2D::cropLineToConvexPolygon(vanishingLineA, viewportAndLocalPoly, true, true);
                KisAlgebra2D::cropLineToConvexPolygon(vanishingLineB, viewportAndLocalPoly, true, true);
                path.moveTo(vanishingLineA.p1());
                path.lineTo(vanishingLineA.p2());
                path.moveTo(vanishingLineB.p1());
                path.lineTo(vanishingLineB.p2());
            }

            if (assistantVisible && !isEditing) {
                drawX(gc, initialTransform.map(p1));
                drawX(gc, initialTransform.map(p2));
                drawX(gc, initialTransform.map(p3));
            }

            // Draw orthocenter and cone of vision
            if (isValid()) {

                // Transform document space-->local persp space
                // Makes calculations easier
                const QTransform t = localTransform(p1,p2,p3);
                const QPointF vp_a = t.map(p1);
                const QPointF vp_b = t.map(p2);

                // Transforms local space-->screen space
                // Used when drawing
                const QTransform inv = t.inverted()*initialTransform;
                const QPointF ortho = orthocenter(vp_a, vp_b, t.map(p3));
                drawX(gc, inv.map(ortho));

                const qreal mid = QLineF(vp_a, vp_b).pointAt(0.5).x();
                const qreal radius = QLineF(vp_a, vp_b).length() / 2.0;
                const qreal sp_distance = sqrt(radius*radius - mid*mid);
                const qreal theta = vp_a.y() - ortho.y();
                const qreal cov_size = sqrt(sp_distance*sp_distance - theta*theta);

                if (isEditing) {

                    // draw 90 degree cone of vision (cov)
                    const qreal actual_size = inv.map(QLineF(ortho, QPointF(cov_size, ortho.y()))).length();
                    const QPointF cov_corner = QPointF(actual_size,actual_size);
                    const QRectF cov_circle = QRectF(inv.map(ortho) - cov_corner, inv.map(ortho) + cov_corner);

                    path.moveTo(inv.map(ortho - QPointF(0,cov_size)));
                    path.lineTo(inv.map(ortho + QPointF(0,cov_size)));
                    path.moveTo(inv.map(ortho - QPointF(cov_size, 0)));
                    path.lineTo(inv.map(ortho + QPointF(cov_size, 0)));
                    path.moveTo(inv.map(ortho));
                    path.addEllipse(cov_circle);

                    // station point guide
                    const QPointF sp = inv.map(QPointF(0,vp_a.y()+sp_distance));
                    path.moveTo(sp);
                    path.lineTo(initialTransform.map(p1));
                    path.moveTo(sp);
                    path.lineTo(initialTransform.map(p2));

                    // vertical tilt guide
                    const QPointF vsp = inv.map(QPointF(ortho + QPointF(cov_size,0)));;
                    const QPointF pp = inv.map(QPointF(0,vp_a.y()));
                    path.moveTo(vsp);
                    path.lineTo(pp);
                    path.moveTo(vsp);
                    path.lineTo(initialTransform.map(p3));
                }

                drawPath(gc, path, isSnappingActive(), true);
                path = QPainterPath();

                // Draw gridlines
                const QPointF principal_pt = QPointF(0,vp_a.y());
                const QPointF cov_edge = QPointF(cov_size, 0);


                Q_FOREACH(QPointF vertical_sp, // Use both left and right vertical SP
                          QList<QPointF>({ortho + cov_edge, ortho - cov_edge})) {

                    // Represents line to the edge of 90 degree cone of vision
                    QLineF line = QLineF(vertical_sp, principal_pt);
                    line.setAngle(line.angle() + 45);

                    // Inner angle between vertical_sp--->ortho and vertical_sp--->edge of cov, as a value between 0 and 180
                    const qreal vert_angle = qRadiansToDegrees(qAcos(qCos(qDegreesToRadians(line.angleTo(QLineF(vertical_sp, ortho))))));

                    // Represents point on floor/ceiling.
                    // It appears visually fixed when turning left/right
                    QPointF base_pt = QPointF();
                    line.intersect(QLineF(principal_pt, ortho), &base_pt);
                    drawX(gc, inv.map(base_pt));

                    qreal cov_edge_dst = sqrt(2*cov_size*cov_size); // sqrt(a^2 + b^2)
                    qreal base_pt_dst = cov_edge_dst * qCos(qDegreesToRadians(vert_angle));
                    qreal visual_angle = qRadiansToDegrees(qAtan( cov_size / base_pt_dst ));
                    const qreal cone_visual_radius = cov_size * qTan(qDegreesToRadians(visual_angle));

                    const QPointF sp = QPointF(0,vp_a.y()+sp_distance);
                    using VPOpp = QList<QPointF>;
                    Q_FOREACH(VPOpp vps, QList<VPOpp>({{vp_a,vp_b},{vp_b,vp_a}})) {
                        const QPointF horizontal_vp = vps[0];
                        const QPointF horizontal_vp_opp = vps[1];
                        const qreal horz_angle = qRadiansToDegrees(qAsin(qSin(qDegreesToRadians(QLineF(sp,horizontal_vp).angle()))));
                        const qreal interval = cone_visual_radius / qCos(qDegreesToRadians(horz_angle));
                        const QPointF interval_vec = QPointF(interval, 0);
                        path.moveTo(inv.map(base_pt));
                        path.lineTo(inv.map(horizontal_vp));
                        for (int i = -10; i <= 10; i++) {
                            const QPointF next_base_pt = base_pt + i * interval_vec;
                            QLineF gridline = inv.map(QLineF(next_base_pt, horizontal_vp_opp));
                            KisAlgebra2D::cropLineToConvexPolygon(gridline, viewportAndLocalPoly, true, true);
                            if (vert_angle < 90) {
                                path.moveTo(gridline.p1());
                                path.lineTo(inv.map(horizontal_vp_opp));
                            } else {
                                path.moveTo(gridline.p2());
                                path.lineTo(inv.map(horizontal_vp_opp));
                            }

                        }

                    }

                    drawX(gc, inv.map(vertical_sp));
                }

                drawPath(gc, path, isSnappingActive(), true);
            } else {
                drawError(gc, path);
            }
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
    }

    const bool no_intersection =
        // NB: opp_snap_point is initialized here
        sp_to_opp_vp.intersect(horizon, &opp_snap_point) == QLineF::NoIntersection;

    *dragged_handle = inv.map(snap_point);
    *handleOpp = inv.map(opp_snap_point);
    *adjustment = *dragged_handle; // clear
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
