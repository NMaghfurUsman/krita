/*
 * SPDX-FileCopyrightText: 2008 Cyrille Berger <cberger@cberger.net>
 * SPDX-FileCopyrightText: 2010 Geoffry Song <goffrie@gmail.com>
 * SPDX-FileCopyrightText: 2021 Nabil Maghfur Usman <nmaghfurusman@gmail.com>
 *
 *  SPDX-License-Identifier: LGPL-2.0-or-later
 */

#ifndef _THREE_POINT_ASSISTANT_H_
#define _THREE_POINT_ASSISTANT_H_

#include <QObject>
#include <QPolygonF>
#include <QLineF>
#include <QTransform>
#include "kis_painting_assistant.h"
#include "NPointPerspective.h"

class ThreePointAssistant : public NPointPerspective
{
public:
    enum ThreePointHandle {
        FirstHandle,
        SecondHandle,
        VerticalHandle,
        LocalFirstHandle,
        LocalSecondHandle
    };


    ThreePointAssistant();
    QPointF adjustPosition(const QPointF& point, const QPointF& strokeBegin, bool snapToAny) override;
    void adjustLine(QPointF &point, QPointF& strokeBegin) override;

    void endStroke() override;
    KisPaintingAssistantSP clone(QMap<KisPaintingAssistantHandleSP, KisPaintingAssistantHandleSP> &handleMap) const override;

    QPointF getDefaultEditorPosition() const override;
    int numHandles() const override { return isLocal() ? 5 : 3; }

    bool isAssistantComplete() const override;
    bool canBeLocal() const override;

    void realignSideHandles(KisPaintingAssistantHandleSP dragged_handle) override;
    void realignVanishingPoint(KisPaintingAssistantHandleSP dragged_handle, KoPointerEvent* event, QPointF* drag_start, QPointF* adjustment) override;
    void initSideHandles() override;

    QTransform localTransform(QPointF vp_a, QPointF vp_b, QPointF pt_ct);

protected:
    void drawAssistant(QPainter& gc, const QRectF& updateRect, const KisCoordinatesConverter* converter, bool  cached = true,KisCanvas2* canvas=nullptr, bool assistantVisible=true, bool previewVisible=true) override;
    void drawCache(QPainter& gc, const KisCoordinatesConverter *converter,  bool assistantVisible=true) override;

    KisPaintingAssistantHandleSP firstLocalHandle() const override;
    KisPaintingAssistantHandleSP secondLocalHandle() const override;


private:
    explicit ThreePointAssistant(const ThreePointAssistant &rhs, QMap<KisPaintingAssistantHandleSP, KisPaintingAssistantHandleSP> &handleMap);

    QPointF project(const QPointF& pt, const QPointF& strokeBegin, bool snapToAny);
    KisCanvas2 *m_canvas {nullptr};

    bool isValid();

    QLineF m_snapLine;
    double m_gridDensity {1.0};
    bool m_useVertical {true};

    int m_lastUsedPoint {-1}; // last used vanishing point

};

class ThreePointAssistantFactory : public KisPaintingAssistantFactory
{
public:
    ThreePointAssistantFactory();
    ~ThreePointAssistantFactory() override;
    QString id() const override;
    QString name() const override;
    KisPaintingAssistant* createPaintingAssistant() const override;
};

#endif
