/*
 *  SPDX-FileCopyrightText: 2007 Sven Langkamp <sven.langkamp@gmail.com>
 *
 *  SPDX-License-Identifier: GPL-2.0-or-later
 */

#include "kis_shape_selection_model.h"
#include "kis_debug.h"

#include <KoShapeContainer.h>
#include <KoShapeBackground.h>
#include <KoShapeManager.h>

#include "kis_shape_selection.h"
#include "kis_selection.h"
#include "kis_image.h"
#include "kis_update_selection_job.h"


KisShapeSelectionModel::KisShapeSelectionModel(KisImageResolutionProxySP resolutionProxy, KisSelectionWSP selection, KisShapeSelection* shapeSelection)
    : m_resolutionProxy(resolutionProxy)
    , m_parentSelection(selection)
    , m_shapeSelection(shapeSelection)
    , m_updatesEnabled(true)
{
}

KisShapeSelectionModel::~KisShapeSelectionModel()
{
    m_parentSelection = 0;
}

void KisShapeSelectionModel::requestUpdate(const QRect &updateRect)
{
    m_shapeSelection->recalculateOutlineCache();

    if (m_updatesEnabled) {
        m_parentSelection->requestCompressedProjectionUpdate(updateRect);
    }
}

void KisShapeSelectionModel::setResolutionProxy(KisImageResolutionProxySP newResolutionProxy)
{
    const bool resolutionChanged = !m_resolutionProxy->compareResolution(*newResolutionProxy);

    m_resolutionProxy = newResolutionProxy;

    if (resolutionChanged) {
        requestUpdate(QRect());
    }
}

KisImageResolutionProxySP KisShapeSelectionModel::resolutionProxy() const
{
    return m_resolutionProxy;
}

void KisShapeSelectionModel::add(KoShape *child)
{
    if (!m_shapeSelection) return;

    if (m_shapeMap.contains(child))
        return;

    child->setStroke(KoShapeStrokeModelSP());
    child->setBackground( QSharedPointer<KoShapeBackground>(0));
    m_shapeMap.insert(child, child->boundingRect());
    m_shapeSelection->shapeManager()->addShape(child);

    QRect updateRect = child->boundingRect().toAlignedRect();
    QTransform matrix;
    matrix.scale(m_resolutionProxy->xRes(), m_resolutionProxy->yRes());
    updateRect = matrix.mapRect(updateRect);

    if (m_shapeMap.count() == 1) {
        // The shape is the first one, so the shape selection just got created
        // Pixel selection provides no longer the datamanager of the selection
        // so update the whole selection
        requestUpdate(QRect());
    } else {
        requestUpdate(updateRect);
    }
}

void KisShapeSelectionModel::remove(KoShape *child)
{
    if (!m_shapeMap.contains(child)) return;

    QRect updateRect = child->boundingRect().toAlignedRect();
    m_shapeMap.remove(child);

    if (m_shapeSelection) {
        m_shapeSelection->shapeManager()->remove(child);
    }
    QTransform matrix;
    matrix.scale(m_resolutionProxy->xRes(), m_resolutionProxy->yRes());
    updateRect = matrix.mapRect(updateRect);
    if (m_shapeSelection) { // No m_shapeSelection indicates the selection is being deleted
        requestUpdate(updateRect);
    }
}

void KisShapeSelectionModel::setUpdatesEnabled(bool enabled)
{
    m_updatesEnabled = enabled;
}

bool KisShapeSelectionModel::updatesEnabled() const
{
    return m_updatesEnabled;
}

void KisShapeSelectionModel::setClipped(const KoShape *child, bool clipping)
{
    Q_UNUSED(child);
    Q_UNUSED(clipping);
}

bool KisShapeSelectionModel::isClipped(const KoShape *child) const
{
    Q_UNUSED(child);
    return false;
}

void KisShapeSelectionModel::setInheritsTransform(const KoShape *shape, bool inherit)
{
    Q_UNUSED(shape);
    Q_UNUSED(inherit);
}

bool KisShapeSelectionModel::inheritsTransform(const KoShape *shape) const
{
    Q_UNUSED(shape);
    return false;
}

int KisShapeSelectionModel::count() const
{
    return m_shapeMap.count();
}

QList<KoShape*> KisShapeSelectionModel::shapes() const
{
    return QList<KoShape*>(m_shapeMap.keys());
}
void KisShapeSelectionModel::containerChanged(KoShapeContainer *, KoShape::ChangeType)
{
}

void KisShapeSelectionModel::childChanged(KoShape * child, KoShape::ChangeType type)
{
    if (!m_shapeSelection) return;

    // TODO: check if still needed
    if (type == KoShape::ParentChanged) return;

    QRectF changedRect = m_shapeMap[child];
    changedRect = changedRect.united(child->boundingRect());
    m_shapeMap[child] = child->boundingRect();

    QTransform matrix;
    matrix.scale(m_resolutionProxy->xRes(), m_resolutionProxy->yRes());
    changedRect = matrix.mapRect(changedRect);

    requestUpdate(changedRect.toAlignedRect());
}

void KisShapeSelectionModel::setShapeSelection(KisShapeSelection* selection)
{
    m_shapeSelection = selection;
}
