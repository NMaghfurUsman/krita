/*
 * SPDX-FileCopyrightText: 2008 Cyrille Berger <cberger@cberger.net>
 * SPDX-FileCopyrightText: 2010 Geoffry Song <goffrie@gmail.com>
 * SPDX-FileCopyrightText: 2021 Nabil Maghfur Usman <nmaghfurusman@gmail.com>
 *
 *  SPDX-License-Identifier: LGPL-2.0-or-later
 */

#include "TwoPointAssistant.h"
#include "kis_debug.h"
#include "kis_types.h"
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

TwoPointAssistant::TwoPointAssistant()
    : NPointPerspective("two point", i18n("Two point assistant"))
{
}

TwoPointAssistant::TwoPointAssistant(const TwoPointAssistant &rhs, QMap<KisPaintingAssistantHandleSP, KisPaintingAssistantHandleSP> &handleMap)
    : NPointPerspective(rhs, handleMap)
    , m_canvas(rhs.m_canvas)
    , m_snapLine(rhs.m_snapLine)
    , m_gridDensity(rhs.m_gridDensity)
    , m_useVertical(rhs.m_useVertical)
    , m_lastUsedPoint(rhs.m_lastUsedPoint)
{
}

KisPaintingAssistantSP TwoPointAssistant::clone(QMap<KisPaintingAssistantHandleSP, KisPaintingAssistantHandleSP> &handleMap) const
{
    return KisPaintingAssistantSP(new TwoPointAssistant(*this, handleMap));
}

QPointF TwoPointAssistant::project(const QPointF& point, const QPointF& strokeBegin, const bool snapToAny)
{
    Q_ASSERT(isAssistantComplete());

    QPointF best_pt = point;
    double best_dist = DBL_MAX;
    QList<int> possibleHandles;

    bool overrideLastUsedPoint = false;
    bool useBeginInstead = false;

    // must be above or equal to 0;
    // if useVertical, then last used point must be below 3, because 2 means vertical
    //     and it's the last possible point here (sanity check)
    // if !useVertical, then it must be below 2, because 2 means vertical
    bool isLastUsedPointCorrectNow = m_lastUsedPoint >= 0 && (m_useVertical ? m_lastUsedPoint < 3 : m_lastUsedPoint < 2);

    if (isLocal() && handles().size() == 5) {
        // here we can just return since we don't want to do anything
        // so we're returning a NaN
        // but only if we don't have a point/axes it was already using

        QRectF rect = getLocalRect();
        bool insideLocalRect = rect.contains(point);
        if (!insideLocalRect && (!isLastUsedPointCorrectNow || !m_hasBeenInsideLocalRect)) {
            return QPointF(qQNaN(), qQNaN());
        } else if (insideLocalRect) {
            m_hasBeenInsideLocalRect = true;
        }
    }


    if (!snapToAny && isLastUsedPointCorrectNow) {
        possibleHandles = QList<int>({m_lastUsedPoint});
    } else {
        if (m_useVertical) {
            possibleHandles = QList<int>({0, 1, 2});
        } else {
            possibleHandles = QList<int>({0, 1});
        }
        overrideLastUsedPoint = true;
    }

    Q_FOREACH (int vpIndex, possibleHandles) {
        QPointF vp = *handles()[vpIndex];
        double dist = 0;
        QPointF pt = QPointF();
        QLineF snapLine = QLineF();

        // TODO: Would be a good idea to generalize this whole routine
        // in KisAlgebra2d, as it's all lifted from the vanishing
        // point assistant and parallel ruler assistant, and by
        // extension the perspective assistant...
        qreal dx = point.x() - strokeBegin.x();
        qreal dy = point.y() - strokeBegin.y();

        if (dx * dx + dy * dy < 4.0) {
            // we cannot return here because m_lastUsedPoint needs to be set properly
            useBeginInstead = true;
        }

        if (vp != *handles()[2]) {
            snapLine = QLineF(vp, strokeBegin);
        } else {
            QLineF vertical = QLineF(*handles()[0],*handles()[1]).normalVector();
            snapLine = QLineF(vertical.p1(), vertical.p2());
            QPointF translation = (vertical.p1()-strokeBegin)*-1.0;
            snapLine = snapLine.translated(translation);
        }

        dx = snapLine.dx();
        dy = snapLine.dy();

        const qreal dx2 = dx * dx;
        const qreal dy2 = dy * dy;
        const qreal invsqrlen = 1.0 / (dx2 + dy2);

        pt = QPointF(dx2 * point.x() + dy2 * snapLine.x1() + dx * dy * (point.y() - snapLine.y1()),
                     dx2 * snapLine.y1() + dy2 * point.y() + dx * dy * (point.x() - snapLine.x1()));

        pt *= invsqrlen;
        dist = qAbs(pt.x() - point.x()) + qAbs(pt.y() - point.y());

        if (dist < best_dist) {
            best_pt = pt;
            best_dist = dist;
            if (overrideLastUsedPoint) {
                m_lastUsedPoint = vpIndex;
            }
        }
    }

    return useBeginInstead ? strokeBegin : best_pt;
}

void TwoPointAssistant::endStroke()
{
    m_snapLine = QLineF();
    m_lastUsedPoint = -1;
    KisPaintingAssistant::endStroke();
}

QPointF TwoPointAssistant::adjustPosition(const QPointF& pt, const QPointF& strokeBegin, const bool snapToAny, qreal /*moveThresholdPt*/)
{
    return project(pt, strokeBegin, snapToAny);
}

void TwoPointAssistant::adjustLine(QPointF &point, QPointF &strokeBegin)
{
    QPointF p = project(point, strokeBegin, true);
    point = p;
}

void TwoPointAssistant::drawAssistant(QPainter& gc, const QRectF& updateRect, const KisCoordinatesConverter* converter, bool cached, KisCanvas2* canvas, bool assistantVisible, bool previewVisible)
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

            // Draw circle to represent center of vision
            if (handles().length() == 3 && handle == handles()[2]) {
                const QLineF horizon = QLineF(*handles()[0],*handles()[1]);
                QLineF normal = horizon.normalVector();
                normal.translate(*handles()[2]-normal.p1());
                QPointF cov = horizon.center();
                normal.intersect(horizon,&cov);
                const QPointF center = initialTransform.map(cov);
                QRectF center_ellipse = QRectF(QPointF(center.x() -15, center.y() -15), QSizeF(30, 30));
                QPainterPath pathCenter;
                pathCenter.addEllipse(center_ellipse);
                drawPath(gc, pathCenter, isSnappingActive());
            }
        }

        if (handles().size() <= 2) {
            QPainterPath path;
            int tempDensity = m_gridDensity * 10; // the vanishing point density seems visibly more dense, hence let's make it less dense
            QRect viewport = gc.viewport();

            for (int i = 0; i < handles().size(); i++) {
                const QPointF p = initialTransform.map(*handles()[i]);
                for (int currentAngle=0; currentAngle <= 180; currentAngle = currentAngle + tempDensity) {

                    // determine the correct angle based on the iteration
                    float xPos = cos(currentAngle * M_PI / 180);
                    float yPos = sin(currentAngle * M_PI / 180);
                    float length = 100;
                    QPointF unit = QPointF(length*xPos, length*yPos);

                    // find point
                    QLineF snapLine = QLineF(p, p + unit);
                    if (KisAlgebra2D::intersectLineRect(snapLine, viewport, false)) {
                        // make a line from VP center to edge of canvas with that angle

                        path.moveTo(snapLine.p1());
                        path.lineTo(snapLine.p2());
                    }

                    QLineF snapLine2 = QLineF(p, p - unit);
                    if (KisAlgebra2D::intersectLineRect(snapLine2, viewport, false)) {
                        // make a line from VP center to edge of canvas with that angle

                        path.moveTo(snapLine2.p1());
                        path.lineTo(snapLine2.p2());
                    }


                }

                drawPreview(gc, path);//and we draw the preview.

            }
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

        // draw the side handle bars
        if (isEditing == true && !sideHandles().isEmpty()) {
            path.moveTo(initialTransform.map(p1));
            path.lineTo(initialTransform.map(*sideHandles()[0]));
            path.lineTo(initialTransform.map(*sideHandles()[1]));
            path.moveTo(initialTransform.map(p2));
            path.lineTo(initialTransform.map(*sideHandles()[2]));
            path.lineTo(initialTransform.map(*sideHandles()[3]));
            path.moveTo(initialTransform.map(p1));
            path.lineTo(initialTransform.map(*sideHandles()[4]));
            path.lineTo(initialTransform.map(*sideHandles()[5]));
            path.moveTo(initialTransform.map(p2));
            path.lineTo(initialTransform.map(*sideHandles()[6]));
            path.lineTo(initialTransform.map(*sideHandles()[7]));
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


        drawPreview(gc,previewPath);
        drawPath(gc, path, isSnappingActive());

        if (handles().size() >= 3 && isSnappingActive()) {
            path = QPainterPath(); // clear
            const QPointF p3 = *handles()[2];

            qreal size = 0;
            const QTransform t = localTransform(p1,p2,p3,&size);
            const QTransform inv = t.inverted();
            const QPointF vp_a = t.map(p1);
            const QPointF vp_b = t.map(p2);

            if ((vp_a.x() < 0 && vp_b.x() > 0) ||
                (vp_a.x() > 0 && vp_b.x() < 0)) {
                if (m_useVertical) {
                    // Draw vertical line, but only if the center is between both VPs
                    QLineF vertical = initialTransform.map(inv.map(QLineF::fromPolar(1,90)));
                    if (!isEditing) vertical.translate(mousePos - vertical.p1());
                    KisAlgebra2D::cropLineToConvexPolygon(vertical, viewportAndLocalPoly, true, true);
                    if (previewVisible) {
                        path.moveTo(vertical.p1());
                        path.lineTo(vertical.p2());
                    }

                    if (assistantVisible) {
                        // Display a notch to represent the center of vision
                        path.moveTo(initialTransform.map(inv.map(QPointF(0,vp_a.y()-10))));
                        path.lineTo(initialTransform.map(inv.map(QPointF(0,vp_a.y()+10))));
                    }
                    drawPreview(gc,path);
                    path = QPainterPath(); // clear
                }
            }

            const QPointF upper = QPointF(0,vp_a.y() + size);
            const QPointF lower = QPointF(0,vp_a.y() - size);

            // Set up the fading effect for the grid lines
            // Needed so the grid density doesn't look distracting
            QColor color = effectiveAssistantColor();
            QGradient fade = QLinearGradient(initialTransform.map(inv.map(upper)),
                                             initialTransform.map(inv.map(lower)));
            color.setAlphaF(0);
            fade.setColorAt(0.4, effectiveAssistantColor());
            fade.setColorAt(0.5, color);
            fade.setColorAt(0.6, effectiveAssistantColor());
            const QPen pen = gc.pen();
            const QBrush new_brush = QBrush(fade);
            int width = 1;
            const QPen new_pen = QPen(new_brush, width, pen.style());
            gc.setPen(new_pen);

            const QList<QPointF> station_points = {upper, lower};
            const QList<QPointF> vanishing_points = {vp_a, vp_b};

            // Draw grid lines above and below the horizon
            Q_FOREACH (const QPointF sp, station_points) {

                // Draw grid lines towards each vanishing point
                Q_FOREACH (const QPointF vp, vanishing_points) {

                    // Interval between each grid line, uses grid density specified by user
                    const qreal initial_angle = QLineF(sp, vp).angle();
                    const qreal interval = size*m_gridDensity / cos((initial_angle - 90) * M_PI/180);
                    const QPointF translation = QPointF(interval, 0);

                    // Draw grid lines originating from both the left and right of the central vertical line
                    Q_FOREACH (const int dir, QList<int>({-1, 1})) {

                        // Limit at 300 grid lines per direction, reasonable even for m_gridDensity=0.1;
                        for (int i = 0; i <= 300; i++) {
                            const QLineF gridline = QLineF(sp + translation * i * dir, vp);

                            // Don't bother drawing lines that are nearly parallel to horizon
                            const qreal angle = gridline.angle();
                            if (angle < 0.25 || angle > 359.75 || (angle < 180.25 && angle > 179.75)) {
                                break;
                            }

                            QLineF drawn_gridline = initialTransform.map(inv.map(gridline));
                            KisAlgebra2D::cropLineToConvexPolygon(drawn_gridline, viewportAndLocalPoly, true, false);

                            if (assistantVisible || isEditing == true) {
                                path.moveTo(drawn_gridline.p2());
                                path.lineTo(drawn_gridline.p1());
                            }
                        }
                    }
                }
            }
            gc.drawPath(path);
        }
    }

    gc.restore();
    //KisPaintingAssistant::drawAssistant(gc, updateRect, converter, cached, canvas, assistantVisible, previewVisible);
}

void TwoPointAssistant::drawCache(QPainter& gc, const KisCoordinatesConverter *converter, bool assistantVisible)
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

KisPaintingAssistantHandleSP TwoPointAssistant::firstLocalHandle() const
{
    if (handles().size() > LocalFirstHandle) {
        return handles().at(LocalFirstHandle);
    } else {
        return nullptr;
    }
}

KisPaintingAssistantHandleSP TwoPointAssistant::secondLocalHandle() const
{
    if (handles().size() > LocalSecondHandle) {
        return handles().at(LocalSecondHandle);
    } else {
        return nullptr;
    }
}

QPointF TwoPointAssistant::getDefaultEditorPosition() const
{
    int centerOfVisionHandle = 2;
    if (handles().size() > centerOfVisionHandle) {
        return *handles().at(centerOfVisionHandle);
    } else if (handles().size() > 0) {
        KIS_SAFE_ASSERT_RECOVER_RETURN_VALUE(false, *handles().at(0));
        return *handles().at(0);
    } else {
        KIS_SAFE_ASSERT_RECOVER_RETURN_VALUE(false, QPointF(0, 0));
        return QPointF(0, 0);
    }
}

void TwoPointAssistant::setGridDensity(double density)
{
    m_gridDensity = density;
}

bool TwoPointAssistant::useVertical()
{
    return m_useVertical;
}

void TwoPointAssistant::setUseVertical(bool value)
{
    m_useVertical = value;
}

double TwoPointAssistant::gridDensity()
{
    return m_gridDensity;
}

QTransform TwoPointAssistant::localTransform(QPointF vp_a, QPointF vp_b, QPointF pt_c, qreal* size)
{
    QTransform t = QTransform();
    t.rotate(QLineF(vp_a, vp_b).angle());
    t.translate(-pt_c.x(),-pt_c.y());
    const QLineF horizon = QLineF(t.map(vp_a), QPointF(t.map(vp_b).x(),t.map(vp_a).y()));
    *size = sqrt(pow(horizon.length()/2.0,2) - pow(abs(horizon.center().x()),2));

    return t;
}

bool TwoPointAssistant::isAssistantComplete() const
{
    return handles().size() >= numHandles();
}

bool TwoPointAssistant::canBeLocal() const
{
    return true;
}

void TwoPointAssistant::saveCustomXml(QXmlStreamWriter* xml)
{
    xml->writeStartElement("gridDensity");
    xml->writeAttribute("value", KisDomUtils::toString( this->gridDensity()));
    xml->writeEndElement();
    xml->writeStartElement("useVertical");
    xml->writeAttribute("value", KisDomUtils::toString( (int)this->useVertical()));
    xml->writeEndElement();
    xml->writeStartElement("isLocal");
    xml->writeAttribute("value", KisDomUtils::toString( (int)this->isLocal()));
    xml->writeEndElement();

}

bool TwoPointAssistant::loadCustomXml(QXmlStreamReader* xml)
{
    if (xml && xml->name() == "gridDensity") {
        this->setGridDensity((float)KisDomUtils::toDouble(xml->attributes().value("value").toString()));
    }
    if (xml && xml->name() == "useVertical") {
        this->setUseVertical((bool)KisDomUtils::toInt(xml->attributes().value("value").toString()));
    }
    if (xml && xml->name() == "isLocal") {
        this->setLocal((bool)KisDomUtils::toInt(xml->attributes().value("value").toString()));
    }
    return true;
}

void TwoPointAssistant::realignSideHandles(KisPaintingAssistantHandleSP dragged_handle) {

    const bool far_handle_is_dragged =
        dragged_handle == sideHandles()[1] ||
        dragged_handle == sideHandles()[3] ||
        dragged_handle == sideHandles()[5] ||
        dragged_handle == sideHandles()[7];

    if (far_handle_is_dragged) {
        QLineF perspective_line_a, perspective_line_b;
        QPointF vp_new_pos(0, 0);
        KisPaintingAssistantHandleSP vp_moved;
        if (dragged_handle == sideHandles()[1] ||
            dragged_handle == sideHandles()[5]) {
            vp_moved = handles()[0];
            perspective_line_a = QLineF(*sideHandles()[0], *sideHandles()[1]);
            perspective_line_b = QLineF(*sideHandles()[4], *sideHandles()[5]);
        } else {
            vp_moved = handles()[1];
            perspective_line_a = QLineF(*sideHandles()[3], *sideHandles()[2]);
            perspective_line_b = QLineF(*sideHandles()[6], *sideHandles()[7]);
        }
        if (perspective_line_a.intersect(perspective_line_b,
                                         &vp_new_pos) !=
            QLineF::NoIntersection) {
            *vp_moved = vp_new_pos;
        }
    } else {
        QLineF perspective_line_a1;
        QLineF perspective_line_b1;
        QLineF perspective_line_a2;
        QLineF perspective_line_b2;

        perspective_line_a1 = QLineF(*handles()[0], *sideHandles()[0]);
        perspective_line_a1.setLength(
            QLineF(*sideHandles()[0], *sideHandles()[1]).length());
        perspective_line_a1.translate(*sideHandles()[0] -
                                      perspective_line_a1.p1());
        *sideHandles()[1] = perspective_line_a1.p2();

        perspective_line_b1 = QLineF(*handles()[0], *sideHandles()[4]);
        perspective_line_b1.setLength(
            QLineF(*sideHandles()[4], *sideHandles()[5]).length());
        perspective_line_b1.translate(*sideHandles()[4] -
                                      perspective_line_b1.p1());
        *sideHandles()[5] = perspective_line_b1.p2();

        perspective_line_a2 = QLineF(*handles()[1], *sideHandles()[2]);
        perspective_line_a2.setLength(
            QLineF(*sideHandles()[2], *sideHandles()[3]).length());
        perspective_line_a2.translate(*sideHandles()[2] -
                                      perspective_line_a2.p1());
        *sideHandles()[3] = perspective_line_a2.p2();

        perspective_line_b2 = QLineF(*handles()[1], *sideHandles()[6]);
        perspective_line_b2.setLength(
            QLineF(*sideHandles()[6], *sideHandles()[7]).length());
        perspective_line_b2.translate(*sideHandles()[6] -
                                      perspective_line_b2.p1());
        *sideHandles()[7] = perspective_line_b2.p2();
    }
}

void TwoPointAssistant::realignVanishingPoint(KisPaintingAssistantHandleSP dragged_handle, KoPointerEvent* event, QPointF* drag_start, QPointF* adjustment)
{
    // Snapping interactions that are specific to the two point assistant.
    // Skip this code block when only Shift is pressed, as
    // Shift means we only need closest-axis snapping.
    KisPaintingAssistantHandleSP handleOpp = dragged_handle == handles()[0] ? handles()[1] : handles()[0];
    const QPointF prevPoint = adjustment->isNull() ? *drag_start : *adjustment;

    qreal size = 0;
    const QTransform t = localTransform(prevPoint,*handleOpp,*handles()[2],&size);
    const QTransform inv = t.inverted();

    // Exact alignment matters here, so fudge horizon line
    // to be perfectly horizontal instead of trusting the
    // QTransform calculation to do it
    const QLineF horizon = QLineF(t.map(prevPoint), QPointF(t.map(*handleOpp).x(),t.map(prevPoint).y()));
    const QPointF sp = QPointF(0,horizon.p1().y()+size);

    const bool preserve_distortion_snap = event->modifiers() == Qt::ControlModifier;
    const bool preserve_left_right_ratio_snap = event->modifiers() == (Qt::ControlModifier|Qt::ShiftModifier);
    const bool preserve_horizon_snap = event->modifiers() == Qt::AltModifier;

    QPointF snap_point;
    QPointF opp_snap_point;
    QLineF sp_to_opp_vp;

    if (preserve_distortion_snap) {
        const QLineF sp_to_vp = QLineF(sp, t.map(*dragged_handle));
        sp_to_opp_vp = sp_to_vp.normalVector();
        sp_to_vp.intersect(horizon,&snap_point);
    } else if (preserve_left_right_ratio_snap) {
        const QLineF prev_sp_to_vp = QLineF(sp, horizon.p1());
        QLineF new_sp_to_vp = prev_sp_to_vp.translated(t.map(*dragged_handle)-sp);
        QPointF new_sp;
        new_sp_to_vp.intersect(QLineF(QPoint(0,0),QPointF(0,1)),&new_sp);
        sp_to_opp_vp = new_sp_to_vp.normalVector().translated(new_sp-new_sp_to_vp.p1());
        new_sp_to_vp.intersect(horizon,&snap_point);
    } else if (preserve_horizon_snap) {
        snap_point = QPointF(t.map(*dragged_handle).x(),horizon.p1().y());
        sp_to_opp_vp = QLineF(sp,QPointF(t.map(prevPoint).x(),horizon.p1().y())).normalVector();
    }

    // The snapping modes must be robust against falling into
    // invalid configurations, so test if the new snap points
    // actually do make sense
    const bool no_intersection =
        // NB: opp_snap_point is initialized here
        sp_to_opp_vp.intersect(horizon, &opp_snap_point) == QLineF::NoIntersection;
    const bool origin_is_between =
        (snap_point.x() < 0 && opp_snap_point.x() > 0) ||
        (snap_point.x() > 0 && opp_snap_point.x() < 0);
    const bool null_opp_point =
        qFuzzyIsNull(opp_snap_point.x()) ||
        qFuzzyIsNull(opp_snap_point.y());
    const bool overlapping_snap_points =
        qFuzzyCompare(opp_snap_point.x(),snap_point.x());

    // Revert to original state if new points are invalid
    if (!origin_is_between || no_intersection || null_opp_point || overlapping_snap_points) {
        *dragged_handle = *drag_start;
        QPointF oppStart;
        // Use different recovery method for different
        // snapping modes
        if (preserve_distortion_snap) {
            sp_to_opp_vp = QLineF(sp, t.map(*drag_start)).normalVector();
            sp_to_opp_vp.intersect(horizon, &oppStart);
        } else {
            const QPointF p1 = t.map(*drag_start);
            const qreal p2x = preserve_horizon_snap ? t.map(*handleOpp).x() : -p1.x();
            const QPointF p2 = QPointF(p2x,p1.y());
            const QLineF new_horizon = QLineF(p1,p2);
            const qreal new_size = sqrt(pow(new_horizon.length()/2.0,2) -
                                        pow(abs(new_horizon.center().x()),2));
            const QPointF new_sp = QPointF(0,horizon.p1().y()+new_size);
            sp_to_opp_vp = QLineF(new_sp, t.map(*drag_start)).normalVector();
        }
        sp_to_opp_vp.intersect(horizon, &oppStart);
        *handleOpp=inv.map(oppStart);
        *adjustment = QPointF(0,0); // clear
    } else {
        // otherwise use the new configuration if it's valid
        *dragged_handle = inv.map(snap_point);
        *handleOpp = inv.map(opp_snap_point);
        *adjustment = *dragged_handle; // clear
    }
}

void TwoPointAssistant::initSideHandles()
{
    if (*handles()[0] == *handles()[1] || *handles()[1] == *handles()[2]) {
        // Place handles() in sensible default position if any of
        // them are overlapping (maybe because user
        // double-clicked)
        const QTransform transform = m_canvas->coordinatesConverter()->documentToWidgetTransform();
        const QTransform inverted = transform.inverted();
        const qreal size = inverted.map(QPointF(m_canvas->canvasWidget()->width(),0)).x();
        *handles()[0] = *handles()[2] - QPointF(-size/3.0,0);
        *handles()[1] = *handles()[2] - QPointF(size/3.0,0);
    }

    const QPointF vp1 = *handles()[0];
    const QPointF vp2 = *handles()[1];
    const QPointF middle = *handles()[2];

    qreal size = 0;
    QTransform t = localTransform(vp1,vp2,middle,&size);
    QTransform inv = t.inverted();

    if (t.map(vp1).x() * t.map(vp2).x() > 0) {
        // Put third handle between first and second if user
        // placed it outside of them, then re-define the transform
        const QLineF horizon = QLineF(t.map(vp1),t.map(vp2));
        const QPointF origin = QPointF(horizon.center().x(),0);
        *handles()[2] = inv.map(origin);
        t = localTransform(vp1,vp2,*handles()[2],&size);
        inv = t.inverted();
    }

    const QPointF above = inv.map(QPointF(0,t.map(vp1).y()+size));
    const QPointF below = inv.map(QPointF(0,t.map(vp1).y()-size));

    Q_FOREACH (QPointF side, QList<QPointF>({above,below})) {
        Q_FOREACH (QPointF vp, QList<QPointF>({vp1, vp2})) {
            QLineF bar = QLineF(side, vp);
            addHandle(new KisPaintingAssistantHandle(bar.pointAt(0.8)), HandleType::SIDE);
            addHandle(new KisPaintingAssistantHandle(bar.pointAt(0.4)), HandleType::SIDE);
        }
    }
}

TwoPointAssistantFactory::TwoPointAssistantFactory()
{
}

TwoPointAssistantFactory::~TwoPointAssistantFactory()
{
}

QString TwoPointAssistantFactory::id() const
{
    return "two point";
}

QString TwoPointAssistantFactory::name() const
{
    return i18n("2 Point Perspective");
}

KisPaintingAssistant* TwoPointAssistantFactory::createPaintingAssistant() const
{
    return new TwoPointAssistant;
}
