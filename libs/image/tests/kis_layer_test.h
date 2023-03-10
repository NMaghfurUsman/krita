/*
 *  SPDX-FileCopyrightText: 2007 Boudewijn Rempt <boud@valdyas.org>
 *
 *  SPDX-License-Identifier: GPL-2.0-or-later
 */

#ifndef KIS_LAYER_TESTER_H
#define KIS_LAYER_TESTER_H

#include <simpletest.h>

#include "kis_layer.h"
#include "kis_types.h"
#include "kis_node_visitor.h"
#include "kis_image.h"
#include "kis_default_bounds.h"

class TestLayer : public KisLayer
{

    Q_OBJECT

public:

    TestLayer(KisImageWSP image, const QString & name, quint8 opacity)
            : KisLayer(image, name, opacity) {
        m_original = new KisPaintDevice(this, image->colorSpace(), new KisDefaultBounds(image));
    }

    KisNodeSP clone() {
        return new TestLayer(*this);
    }
    bool allowAsChild(KisNodeSP) const override {
        return true;
    }

    virtual QString nodeType() {
        return "TEST";
    }

    KisPaintDeviceSP original() const override {
        return m_original;
    }

    KisPaintDeviceSP paintDevice() const override {
        return 0;
    }

    QIcon icon() const override {
        return QIcon();
    }

    KisNodeSP clone() const override {
        return new TestLayer(image(), name(), opacity());
    }

    qint32 x() const override {
        return 0;
    }

    void setX(qint32) override {
    }

    qint32 y() const override {
        return 0;
    }

    void setY(qint32) override {
    }

    QRect extent() const override {
        return QRect();
    }

    QRect exactBounds() const override {
        return QRect();
    }

    using KisLayer::accept;

    bool accept(KisNodeVisitor& v) override {
        return v.visit(this);
    }

private:
    KisPaintDeviceSP m_original;
};

class KisLayerTest : public QObject
{
    Q_OBJECT

private Q_SLOTS:

    void testCreation();
    void testOrdering();
    void testMoveNode();
    void testMoveLayer();
    void testMasksChangeRect();
    void testMoveLayerWithMaskThreaded();
};

#endif

