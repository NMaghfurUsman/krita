/*
 *  SPDX-FileCopyrightText: 2010 Adam Celarek <kdedev at xibo dot at>
 *
 *  SPDX-License-Identifier: LGPL-2.0-or-later
 */

#include "kis_color_history.h"
#include "kis_canvas2.h"
#include "KisViewManager.h"
#include "KisView.h"
#include "kis_canvas_resource_provider.h"

#include <KoCompositeOpRegistry.h>

#include <QToolButton>
#include <QColor>
#include <QList>
#include <kis_icon.h>

KisColorHistory::KisColorHistory(QWidget *parent)
    : KisColorPatches("lastUsedColors", parent)
    , m_resourceProvider(0)
{
    m_clearButton = new QToolButton(this);
    m_clearButton->setIcon(KisIconUtils::loadIcon("dialog-cancel-16"));
    m_clearButton->setToolTip(i18n("Clear all color history"));
    m_clearButton->setAutoRaise(true);
    connect(m_clearButton, SIGNAL(clicked()), this, SLOT(clearColorHistory()));

    setAdditionalButtons({m_clearButton});
}

void KisColorHistory::unsetCanvas()
{
    KisColorPatches::unsetCanvas();
    m_resourceProvider = 0;
}

void KisColorHistory::setCanvas(KisCanvas2 *canvas)
{
    if (!canvas) return;

    KisColorPatches::setCanvas(canvas);

    if (m_resourceProvider) {
        m_resourceProvider->disconnect(this);
    }

    m_resourceProvider = canvas->imageView()->resourceProvider();


    connect(canvas->imageView()->resourceProvider(), SIGNAL(sigFGColorUsed(KoColor)),
            this, SLOT(addColorToHistory(KoColor)), Qt::UniqueConnection);
}

KisColorSelectorBase* KisColorHistory::createPopup() const
{
    KisColorHistory* ret = new KisColorHistory();
    ret->setCanvas(m_canvas);
    ret->setColors(colors());
    return ret;
}

void KisColorHistory::addColorToHistory(const KoColor& color)
{
    // don't add color in erase mode. See https://bugs.kde.org/show_bug.cgi?id=298940
    if (m_resourceProvider && m_resourceProvider->currentCompositeOp() == COMPOSITE_ERASE) return;

    KisColorPatches::addColorPatch(color);
}

void KisColorHistory::clearColorHistory() {
    setColors(QList<KoColor>());
}
