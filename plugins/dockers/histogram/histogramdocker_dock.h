/*
 *  SPDX-FileCopyrightText: 2016 Eugene Ingerman <geneing at gmail dot com>
 *
 *  SPDX-License-Identifier: LGPL-2.0-or-later
 */

#ifndef _HISTOGRAM_DOCK_H_
#define _HISTOGRAM_DOCK_H_

#include <QDockWidget>
#include <QPointer>

#include <KoCanvasObserverBase.h>

#include <kis_paint_device.h>
#include <kis_canvas2.h>

class QVBoxLayout;
class KisIdleWatcher;
class KoHistogramProducer;
class HistogramDockerWidget;

class HistogramDockerDock : public QDockWidget, public KoCanvasObserverBase
{
public:
    HistogramDockerDock();

    QString observerName() override { return "HistogramDockerDock"; }
    void setCanvas(KoCanvasBase *canvas) override;
    void unsetCanvas() override;

private:
    QVBoxLayout *m_layout;
    HistogramDockerWidget *m_histogramWidget;
    QPointer<KisCanvas2> m_canvas;
};


#endif
