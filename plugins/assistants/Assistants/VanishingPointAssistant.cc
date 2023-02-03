/*
 * SPDX-FileCopyrightText: 2008 Cyrille Berger <cberger@cberger.net>
 * SPDX-FileCopyrightText: 2010 Geoffry Song <goffrie@gmail.com>
 * SPDX-FileCopyrightText: 2014 Wolthera van Hövell tot Westerflier <griffinvalley@gmail.com>
 * SPDX-FileCopyrightText: 2017 Scott Petrovic <scottpetrovic@gmail.com>
 *
 *  SPDX-License-Identifier: LGPL-2.0-or-later
 */

#include "VanishingPointAssistant.h"

#include "kis_debug.h"
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
#include <qglobal.h>

VanishingPointAssistant::VanishingPointAssistant()
    : NPointPerspective("vanishing point", i18n("Vanishing Point assistant"))
{
}

VanishingPointAssistant::VanishingPointAssistant(const VanishingPointAssistant &rhs, QMap<KisPaintingAssistantHandleSP, KisPaintingAssistantHandleSP> &handleMap)
    : NPointPerspective(rhs, handleMap)
    , m_canvas(rhs.m_canvas)
    , m_referenceLineDensity(rhs.m_referenceLineDensity)
{
}

KisPaintingAssistantSP VanishingPointAssistant::clone(QMap<KisPaintingAssistantHandleSP, KisPaintingAssistantHandleSP> &handleMap) const
{
    return KisPaintingAssistantSP(new VanishingPointAssistant(*this, handleMap));
}

QPointF VanishingPointAssistant::project(const QPointF& pt, const QPointF& strokeBegin)
{
    //Q_ASSERT(handles().size() == 1 || handles().size() == 5);
    //code nicked from the perspective ruler.
    qreal dx = pt.x() - strokeBegin.x();
    qreal dy = pt.y() - strokeBegin.y();

    if (dx * dx + dy * dy < 4.0) {
        // allow some movement before snapping
        return strokeBegin;
    }

    if (isLocal() && isAssistantComplete()) {
        if (getLocalRect().contains(pt)) {
            m_hasBeenInsideLocalRect = true;
        } else if (!m_hasBeenInsideLocalRect) { // isn't inside and wasn't inside before
            return QPointF(qQNaN(), qQNaN());
        }
    }

    //dbgKrita<<strokeBegin<< ", " <<*handles()[0];
    QLineF snapLine = QLineF(*handles()[0], strokeBegin);


    dx = snapLine.dx();
    dy = snapLine.dy();

    const qreal dx2 = dx * dx;
    const qreal dy2 = dy * dy;
    const qreal invsqrlen = 1.0 / (dx2 + dy2);

    QPointF r(dx2 * pt.x() + dy2 * snapLine.x1() + dx * dy * (pt.y() - snapLine.y1()),
              dx2 * snapLine.y1() + dy2 * pt.y() + dx * dy * (pt.x() - snapLine.x1()));

    r *= invsqrlen;
    return r;
}

QPointF VanishingPointAssistant::adjustPosition(const QPointF& pt, const QPointF& strokeBegin, const bool /*snapToAny*/)
{
    return project(pt, strokeBegin);
}

void VanishingPointAssistant::adjustLine(QPointF &point, QPointF &strokeBegin)
{
    point = project(point, strokeBegin);
}

void VanishingPointAssistant::drawAssistant(QPainter& gc, const QRectF& updateRect, const KisCoordinatesConverter* converter, bool cached, KisCanvas2* canvas, bool assistantVisible, bool previewVisible)
{
    // HACK ALERT: side handles aren't saved in old krita versions
    // we need to just add a default position for now if we are loading a vanishing point
    if (sideHandles().isEmpty()) {
        QPointF vpPoint = *handles()[0]; // main vanishing point
        addHandle(new KisPaintingAssistantHandle(vpPoint + QPointF(-70,0)), HandleType::SIDE);
        addHandle(new KisPaintingAssistantHandle(vpPoint + QPointF(-140,0)), HandleType::SIDE);
        addHandle(new KisPaintingAssistantHandle(vpPoint + QPointF(70,0)), HandleType::SIDE);
        addHandle(new KisPaintingAssistantHandle(vpPoint + QPointF(140,0)), HandleType::SIDE);
    }

    gc.save();
    gc.resetTransform();
    QPointF mousePos(0,0);

    if (canvas) {
        //simplest, cheapest way to get the mouse-position//
        mousePos= canvas->canvasWidget()->mapFromGlobal(QCursor::pos());
        m_canvas = canvas;
    }
    else {
        //...of course, you need to have access to a canvas-widget for that.//
        mousePos = QCursor::pos();//this'll give an offset//
        dbgFile<<"canvas does not exist in ruler, you may have passed arguments incorrectly:"<<canvas;
    }

    QRect viewport= gc.viewport();

    QPolygonF viewportAndLocalPoly = (isLocal() && isAssistantComplete()) ?
                QPolygonF(QRectF(viewport)).intersected(converter->documentToWidgetTransform().map(QPolygonF(QRectF(getLocalRect())))) : QPolygonF(QRectF(viewport));

    // draw controls when we are not editing
    if (canvas && canvas->paintingAssistantsDecoration()->isEditingAssistants() == false && isAssistantComplete()) {

        if (isSnappingActive() && previewVisible == true) {
            //don't draw if invalid.

            QTransform initialTransform = converter->documentToWidgetTransform();
            QPointF startPoint = initialTransform.map(*handles()[0]);

            if (m_followBrushPosition && m_adjustedPositionValid) {
                mousePos = initialTransform.map(m_adjustedBrushPosition);
            }

            QLineF snapLine= QLineF(startPoint, mousePos);

            KisAlgebra2D::cropLineToConvexPolygon(snapLine, viewportAndLocalPoly, false, true);

            QPainterPath path;

            path.moveTo(snapLine.p2());
            path.lineTo(snapLine.p1());

            drawPreview(gc, path);//and we draw the preview.

        }
    }




    // editor specific controls display
    if (canvas->paintingAssistantsDecoration()->isEditingAssistants()) {

        // draws a circle around the vanishing point node while editing
        QTransform initialTransform = converter->documentToWidgetTransform();
        QPointF p0 = initialTransform.map(*handles()[0]); // main vanishing point
        QPointF p1 = initialTransform.map(*sideHandles()[0]);
        QPointF p2 = initialTransform.map(*sideHandles()[1]);
        QPointF p3 = initialTransform.map(*sideHandles()[2]);
        QPointF p4 = initialTransform.map(*sideHandles()[3]);


        QRectF ellipse = QRectF(QPointF(p0.x() -15, p0.y() -15), QSizeF(30, 30));

        QPainterPath pathCenter;
        pathCenter.addEllipse(ellipse);
        drawPath(gc, pathCenter, isSnappingActive());

        QColor paintingColor = effectiveAssistantColor();


        // draw the lines connecting the different nodes
        QPen penStyle(paintingColor, 2.0, Qt::SolidLine);

        if (!isSnappingActive()) {
            QColor snappingColor = paintingColor;
            snappingColor.setAlpha(snappingColor.alpha() * 0.2);

            penStyle.setColor(snappingColor);
        }

        gc.save();
        gc.setPen(penStyle);
        gc.drawLine(p0, p1);
        gc.drawLine(p0, p3);
        gc.drawLine(p1, p2);
        gc.drawLine(p3, p4);
        gc.restore();
    }

    QTransform initialTransform = converter->documentToWidgetTransform();

    // draw the local rectangle
    if (assistantVisible && isLocal() && isAssistantComplete()) {
        // limited area rectangle
        QPainterPath path;
        QPointF p1 = *handles()[(int)LocalFirstHandle];
        QPointF p3 = *handles()[(int)LocalSecondHandle];
        QPointF p2 = QPointF(p1.x(), p3.y());
        QPointF p4 = QPointF(p3.x(), p1.y());

        path.moveTo(initialTransform.map(p1));

        path.lineTo(initialTransform.map(p2));
        path.lineTo(initialTransform.map(p3));
        path.lineTo(initialTransform.map(p4));
        path.lineTo(initialTransform.map(p1));
        drawPath(gc, path, isSnappingActive());//and we draw the rectangle
    }


    // draw references guide for vanishing points at specified density
    if (assistantVisible && this->isSnappingActive() ) {

        // cycle through degrees from 0 to 180. We are doing an infinite line, so we don't need to go 360
        QPointF p0 = initialTransform.map(*handles()[0]); // main vanishing point

        for (int currentAngle=0; currentAngle <= 180; currentAngle = currentAngle + m_referenceLineDensity ) {

            // determine the correct angle based on the iteration
            float xPos = cos(currentAngle * M_PI / 180);
            float yPos = sin(currentAngle * M_PI / 180);
            QPointF unitAngle;
            unitAngle.setX(p0.x() + xPos);
            unitAngle.setY(p0.y() + yPos);

            // find point
            QLineF snapLine= QLineF(p0, unitAngle);
            KisAlgebra2D::intersectLineConvexPolygon(snapLine, viewportAndLocalPoly, true, true);

            // make a line from VP center to edge of canvas with that angle
            QPainterPath path;
            path.moveTo(snapLine.p1());
            path.lineTo(snapLine.p2());
            drawPreview(gc, path);//and we draw the preview.
        }
    }


    gc.restore();

    KisPaintingAssistant::drawAssistant(gc, updateRect, converter, cached, canvas, assistantVisible, previewVisible);
}

void VanishingPointAssistant::drawCache(QPainter& gc, const KisCoordinatesConverter *converter, bool assistantVisible)
{
    if (!m_canvas || !isAssistantComplete()) {
        return;
    }

    if (assistantVisible == false ||   m_canvas->paintingAssistantsDecoration()->isEditingAssistants()) {
        return;
    }

    QTransform initialTransform = converter->documentToWidgetTransform();
    QPointF p0 = initialTransform.map(*handles()[0]);

    // draws an "X"
    QPainterPath path;
    path.moveTo(QPointF(p0.x() - 10.0, p0.y() - 10.0));
    path.lineTo(QPointF(p0.x() + 10.0, p0.y() + 10.0));

    path.moveTo(QPointF(p0.x() - 10.0, p0.y() + 10.0));
    path.lineTo(QPointF(p0.x() + 10.0, p0.y() - 10.0));


    drawPath(gc, path, isSnappingActive());
}

KisPaintingAssistantHandleSP VanishingPointAssistant::firstLocalHandle() const
{
    if (handles().size() > LocalFirstHandle) {
        return handles().at(LocalFirstHandle);
    } else {
        return nullptr;
    }
}

KisPaintingAssistantHandleSP VanishingPointAssistant::secondLocalHandle() const
{
    if (handles().size() > LocalSecondHandle) {
        return handles().at(LocalSecondHandle);
    } else {
        return nullptr;
    }
}

QPointF VanishingPointAssistant::getDefaultEditorPosition() const
{
    int pointHandle = 0;
    if (handles().size() > pointHandle) {
        return *handles().at(pointHandle);
    } else {
        KIS_SAFE_ASSERT_RECOVER_RETURN_VALUE(false, QPointF(0, 0));
        return QPointF(0, 0);
    }
}

void VanishingPointAssistant::setReferenceLineDensity(float value)
{
    // cannot have less than 1 degree value
    if (value < 1.0) {
        value = 1.0;
    }

    m_referenceLineDensity = value;
}

float VanishingPointAssistant::referenceLineDensity()
{
    return m_referenceLineDensity;
}

bool VanishingPointAssistant::isAssistantComplete() const
{
    return handles().size() >= numHandles();
}

bool VanishingPointAssistant::canBeLocal() const
{
    return true;
}

void VanishingPointAssistant::saveCustomXml(QXmlStreamWriter* xml)
{
    xml->writeStartElement("angleDensity");
    xml->writeAttribute("value", KisDomUtils::toString( this->referenceLineDensity()));
    xml->writeEndElement();
    xml->writeStartElement("isLocal");
    xml->writeAttribute("value", KisDomUtils::toString( (int)this->isLocal()));
    xml->writeEndElement();
}

bool VanishingPointAssistant::loadCustomXml(QXmlStreamReader* xml)
{
    if (xml && xml->name() == "angleDensity") {
        this->setReferenceLineDensity((float)KisDomUtils::toDouble(xml->attributes().value("value").toString()));
    }
    if (xml && xml->name() == "isLocal") {
        this->setLocal((bool)KisDomUtils::toInt(xml->attributes().value("value").toString()));
    }

    return true;
}

void VanishingPointAssistant::realignSideHandles(KisPaintingAssistantHandleSP dragged_handle) {
    //for inner handles, the outer handle gets translated.
    const bool far_handle_is_dragged =
        dragged_handle == sideHandles()[1] || dragged_handle == sideHandles()[3];

    if (far_handle_is_dragged) {
        QLineF perspective_line_a, perspective_line_b;
        QPointF vp_new_pos(0,0);
        perspective_line_a = QLineF(*sideHandles()[0],*sideHandles()[1]);
        perspective_line_b = QLineF(*sideHandles()[2],*sideHandles()[3]);
        if (perspective_line_a.intersect(perspective_line_b, &vp_new_pos) != QLineF::NoIntersection) {
            *handles()[0] = vp_new_pos;
        }
    } else {
        QLineF perspective_line_a1;
        QLineF perspective_line_b1;

        perspective_line_a1 = QLineF(*handles()[0], *sideHandles()[0]);
        perspective_line_a1.setLength(QLineF(*sideHandles()[0],*sideHandles()[1]).length());
        perspective_line_a1.translate(*sideHandles()[0] - perspective_line_a1.p1());
        *sideHandles()[1] = perspective_line_a1.p2();

        perspective_line_b1 = QLineF(*handles()[0], *sideHandles()[2]);
        perspective_line_b1.setLength(QLineF(*sideHandles()[2],*sideHandles()[3]).length());
        perspective_line_b1.translate(*sideHandles()[2] - perspective_line_b1.p1());
        *sideHandles()[3] = perspective_line_b1.p2();
    }
}

void VanishingPointAssistant::realignVanishingPoint(KisPaintingAssistantHandleSP dragged_handle, KoPointerEvent* event, QPointF* drag_start, QPointF* adjustment) {
    Q_UNIMPLEMENTED();
}

VanishingPointAssistantFactory::VanishingPointAssistantFactory()
{
}

VanishingPointAssistantFactory::~VanishingPointAssistantFactory()
{
}

QString VanishingPointAssistantFactory::id() const
{
    return "vanishing point";
}

QString VanishingPointAssistantFactory::name() const
{
    return i18n("Vanishing Point");
}

KisPaintingAssistant* VanishingPointAssistantFactory::createPaintingAssistant() const
{
    return new VanishingPointAssistant;
}
