/*
 *  SPDX-FileCopyrightText: 2004 Boudewijn Rempt <boud@valdyas.org>
 *  SPDX-FileCopyrightText: 2007 Sven Langkamp <sven.langkamp@gmail.com>
 *
 *
 *  SPDX-License-Identifier: GPL-2.0-or-later
 */

#include "kis_selection.h"

#include "kundo2command.h"

#include "kis_selection_component.h"
#include "kis_pixel_selection.h"
#include "kis_node_graph_listener.h"
#include "kis_node.h"
#include "kis_image.h"

#include "KisImageResolutionProxy.h"
#include "kis_default_bounds.h"
#include "kis_iterator_ng.h"
#include "KisLazyStorage.h"
#include "KisSelectionUpdateCompressor.h"
#include "kis_simple_stroke_strategy.h"
#include "KisDeleteLaterWrapper.h"
#include "kis_command_utils.h"

#include <QReadWriteLock>
#include <QReadLocker>
#include <QWriteLocker>


struct Q_DECL_HIDDEN KisSelection::Private {
    Private(KisSelection *q)
        : isVisible(true),
          shapeSelection(0),
          updateCompressor(q)

    {
    }

    template <typename T>
    static void safeDeleteShapeSelection(T *object, KisSelection *selection);

    // used for forwarding setDirty signals only
    KisNodeWSP parentNode;

    bool isVisible; //false is the selection decoration should not be displayed
    KisImageResolutionProxySP resolutionProxy;
    KisPixelSelectionSP pixelSelection;
    KisSelectionComponent *shapeSelection;
    KisLazyStorage<KisSelectionUpdateCompressor, KisSelection*> updateCompressor;

    /**
     * This lock makes sure that the shape selection is not reincarnated,
     * while some update jobs still access it via KisSelection::updateProjection().
     */
    QReadWriteLock shapeSelectionPointerLock;
};

template <typename T>
void KisSelection::Private::safeDeleteShapeSelection(T *object, KisSelection *selection)
{
    struct ShapeSelectionReleaseStroke : public KisSimpleStrokeStrategy {
        ShapeSelectionReleaseStroke(T *object)
            : KisSimpleStrokeStrategy(QLatin1String("ShapeSelectionReleaseStroke")),
              m_objectWrapper(makeKisDeleteLaterWrapper(object))
        {
            setRequestsOtherStrokesToEnd(false);
            setClearsRedoOnStart(false);
            setNeedsExplicitCancel(true);

            this->enableJob(JOB_FINISH, true, KisStrokeJobData::BARRIER);
            this->enableJob(JOB_CANCEL, true, KisStrokeJobData::BARRIER);
        }

        ~ShapeSelectionReleaseStroke()
        {
            /// it looks like the strategy has not been executed,
            /// the object will leak...
            KIS_SAFE_ASSERT_RECOVER_NOOP(!m_objectWrapper);
        }

        void finishStrokeCallback() override
        {
            m_objectWrapper->deleteLater();
            m_objectWrapper = 0;
        }

        void cancelStrokeCallback() override
        {
            finishStrokeCallback();
        }

    private:
        KisDeleteLaterWrapper<T*> *m_objectWrapper = 0;
    };

    /**
     * Yes, you see it right. We create a fake object, then put it into the GUI
     * events queue using delete later. This object, on destruction creates a
     * stroke, whose only purpose of life is to give birth to another fake
     * object, which will be put into delete-later queue again. This final
     * object, on destruction will finally release the shape selection.
     *
     * If you don't understand that, relax, it is impossible.
     *
     * This trickery is needed to satisfy the following requirements:
     *
     * 1) KisShapeSelection object should be deleted in the GUI thread (all
     * shape manipulations should happen in the GUI thread). Hence we have the
     * last step of wrapping m_shapeSelection into KisDeleteLaterWrapper.
     *
     * 2) KisShapeSelection cannot be deleted before all the updates using this
     * layer/selection has completed its execution. Hence we create of the
     * stroke that uses a BARRIER job of wait until all the updates are
     * finished.
     *
     * 3) We cannot call image->startStroke() from within
     * safeDeleteShapeSelection() directly, because it may be called in the
     * destructor of KisTransactionData, which may theoretically be called from
     * the destructor of KisStrokeStrategy, which will cause a deadlock in the
     * strokes queue. Hence we do the first layer of wrapping into
     * GuiStrokeWrapper.
     */
    struct GuiStrokeWrapper
    {
        GuiStrokeWrapper(KisImageSP image, T *object)
            : m_image(image), m_object(object)
        {
        }

        ~GuiStrokeWrapper()
        {
            KisImageSP image = m_image;

            if (image) {
                KisStrokeId strokeId = image->startStroke(new ShapeSelectionReleaseStroke(m_object));
                image->endStroke(strokeId);
            } else {
                delete m_object;
            }
        }

        KisImageWSP m_image;
        T *m_object;
    };

    if (selection) {
        KisImageSP image = 0;

        KisNodeSP parentNode = selection->parentNode();
        if (parentNode) {
            image = parentNode->image();
        }

        if (image) {
            makeKisDeleteLaterWrapper(new GuiStrokeWrapper(image, object))->deleteLater();
            object = 0;
        }
    }

    if (object) {
        makeKisDeleteLaterWrapper(object)->deleteLater();
        object = 0;
    }
}

struct KisSelection::ChangeShapeSelectionCommand : public KUndo2Command
{
    ChangeShapeSelectionCommand(KisSelectionWSP selection, KisSelectionComponent *shapeSelection)
        : m_selection(selection),
          m_shapeSelection(shapeSelection)
    {
        m_isFlatten = !shapeSelection;
    }

    ~ChangeShapeSelectionCommand() override
    {
        if (m_shapeSelection) {
            Private::safeDeleteShapeSelection(m_shapeSelection, m_selection ? m_selection.data() : 0);
        }

        if (m_reincarnationCommand) {
            Private::safeDeleteShapeSelection(m_reincarnationCommand.take(), m_selection ? m_selection.data() : 0);
        }
    }

    void undo() override
    {
        KIS_SAFE_ASSERT_RECOVER_RETURN(m_selection);

        if (m_reincarnationCommand) {
            m_reincarnationCommand->undo();
        }

        {
            QWriteLocker l(&m_selection->m_d->shapeSelectionPointerLock);
            std::swap(m_selection->m_d->shapeSelection, m_shapeSelection);
        }

        if (!m_isFlatten) {
            m_selection->requestCompressedProjectionUpdate(QRect());
        }
    }

    void redo() override
    {
        KIS_SAFE_ASSERT_RECOVER_RETURN(m_selection);

        if (m_firstRedo) {
            QReadLocker l(&m_selection->m_d->shapeSelectionPointerLock);

            if (bool(m_selection->m_d->shapeSelection) != bool(m_shapeSelection)) {
                m_reincarnationCommand.reset(
                    m_selection->m_d->pixelSelection->reincarnateWithDetachedHistory(m_isFlatten));
            }
            m_firstRedo = false;

        }

        if (m_reincarnationCommand) {
            m_reincarnationCommand->redo();
        }

        {
            QWriteLocker l(&m_selection->m_d->shapeSelectionPointerLock);
            std::swap(m_selection->m_d->shapeSelection, m_shapeSelection);
        }

        if (!m_isFlatten) {
            m_selection->requestCompressedProjectionUpdate(QRect());
        }
    }

private:
    KisSelectionWSP m_selection;
    KisSelectionComponent *m_shapeSelection = 0;
    QScopedPointer<KUndo2Command> m_reincarnationCommand;
    bool m_firstRedo = true;
    bool m_isFlatten = false;
};

KisSelection::KisSelection()
    : KisSelection(nullptr, nullptr)
{
}

KisSelection::KisSelection(KisDefaultBoundsBaseSP defaultBounds, KisImageResolutionProxySP resolutionProxy)
    : m_d(new Private(this))
{
    if (!defaultBounds) {
        defaultBounds = new KisSelectionEmptyBounds(nullptr);
    }

    if (!resolutionProxy) {
        resolutionProxy.reset(new KisImageResolutionProxy(nullptr));
    }

    m_d->resolutionProxy = resolutionProxy;

    m_d->pixelSelection = new KisPixelSelection(defaultBounds, this);
    m_d->pixelSelection->setParentNode(m_d->parentNode);
}

KisSelection::KisSelection(const KisSelection& rhs)
    : KisShared(),
      m_d(new Private(this))
{
    copyFrom(rhs);
}

KisSelection::KisSelection(const KisPaintDeviceSP source, KritaUtils::DeviceCopyMode copyMode,
                           KisDefaultBoundsBaseSP defaultBounds, KisImageResolutionProxySP resolutionProxy)
    : m_d(new Private(this))
{
    if (!defaultBounds) {
        defaultBounds = new KisSelectionEmptyBounds(0);
    }

    m_d->resolutionProxy = resolutionProxy;
    m_d->pixelSelection = new KisPixelSelection(source, copyMode);
    m_d->pixelSelection->setParentSelection(this);
    m_d->pixelSelection->setParentNode(m_d->parentNode);
    m_d->pixelSelection->setDefaultBounds(defaultBounds);
}

KisSelection &KisSelection::operator=(const KisSelection &rhs)
{
    if (&rhs != this) {
        copyFrom(rhs);
    }
    return *this;
}

void KisSelection::copyFrom(const KisSelection &rhs)
{
    m_d->isVisible = rhs.m_d->isVisible;
    m_d->resolutionProxy = rhs.m_d->resolutionProxy;
    m_d->parentNode = 0; // not supposed to be shared

    Q_ASSERT(rhs.m_d->pixelSelection);
    m_d->pixelSelection = new KisPixelSelection(*rhs.m_d->pixelSelection, KritaUtils::CopyAllFrames);
    m_d->pixelSelection->setParentSelection(this);

    QReadLocker l1(&rhs.m_d->shapeSelectionPointerLock);
    QWriteLocker l2(&m_d->shapeSelectionPointerLock);

    if (rhs.m_d->shapeSelection && !rhs.m_d->shapeSelection->isEmpty()) {
        m_d->shapeSelection = rhs.m_d->shapeSelection->clone(this);
        KIS_SAFE_ASSERT_RECOVER_NOOP(m_d->shapeSelection);
        KIS_SAFE_ASSERT_RECOVER(m_d->shapeSelection &&
                                m_d->shapeSelection != rhs.m_d->shapeSelection) {
            m_d->shapeSelection = 0;
        }
    }
    else {
        if (m_d->shapeSelection) {
            Private::safeDeleteShapeSelection(m_d->shapeSelection, this);
            m_d->shapeSelection = 0;
        }
    }
}

KisSelection::~KisSelection()
{
    delete m_d->shapeSelection;
    delete m_d;
}

void KisSelection::setParentNode(KisNodeWSP node)
{
    m_d->parentNode = node;
    m_d->pixelSelection->setParentNode(node);

    // the updates come through the parent image, so all the updates
    // that happened in the meantime are considered "stalled"
    if (node) {
        m_d->updateCompressor->tryProcessStalledUpdate();
    }
}

// for testing purposes only
KisNodeWSP KisSelection::parentNode() const
{
    return m_d->parentNode;
}

bool KisSelection::outlineCacheValid() const
{
    QReadLocker l(&m_d->shapeSelectionPointerLock);
    return m_d->shapeSelection ||
        m_d->pixelSelection->outlineCacheValid();
}

QPainterPath KisSelection::outlineCache() const
{
    QReadLocker l(&m_d->shapeSelectionPointerLock);

    QPainterPath outline;

    if (m_d->shapeSelection) {
        outline += m_d->shapeSelection->outlineCache();
    } else if (m_d->pixelSelection->outlineCacheValid()) {
        outline += m_d->pixelSelection->outlineCache();
    }

    return outline;
}

void KisSelection::recalculateOutlineCache()
{
    QReadLocker l(&m_d->shapeSelectionPointerLock);

    Q_ASSERT(m_d->pixelSelection);

    if (m_d->shapeSelection) {
        m_d->shapeSelection->recalculateOutlineCache();
    } else if (!m_d->pixelSelection->outlineCacheValid()) {
        m_d->pixelSelection->recalculateOutlineCache();
    }
}

bool KisSelection::thumbnailImageValid() const
{
    return m_d->pixelSelection->thumbnailImageValid();
}

void KisSelection::recalculateThumbnailImage(const QColor &maskColor)
{
    m_d->pixelSelection->recalculateThumbnailImage(maskColor);
}

QImage KisSelection::thumbnailImage() const
{
    return m_d->pixelSelection->thumbnailImage();
}

QTransform KisSelection::thumbnailImageTransform() const
{
    return m_d->pixelSelection->thumbnailImageTransform();
}

bool KisSelection::hasNonEmptyPixelSelection() const
{
    return m_d->pixelSelection && !m_d->pixelSelection->isEmpty();
}

bool KisSelection::hasNonEmptyShapeSelection() const
{
    QReadLocker l(&m_d->shapeSelectionPointerLock);
    return m_d->shapeSelection && !m_d->shapeSelection->isEmpty();
}

bool KisSelection::hasShapeSelection() const
{
    QReadLocker l(&m_d->shapeSelectionPointerLock);
    return m_d->shapeSelection;
}

KisPixelSelectionSP KisSelection::pixelSelection() const
{
    return m_d->pixelSelection;
}

KisSelectionComponent* KisSelection::shapeSelection() const
{
    return m_d->shapeSelection;
}

void KisSelection::convertToVectorSelectionNoUndo(KisSelectionComponent* shapeSelection)
{
    KIS_SAFE_ASSERT_RECOVER_RETURN(shapeSelection);

    shapeSelection->setResolutionProxy(m_d->resolutionProxy);
    QScopedPointer<KUndo2Command> cmd(new ChangeShapeSelectionCommand(this, shapeSelection));
    cmd->redo();
}

KUndo2Command *KisSelection::convertToVectorSelection(KisSelectionComponent *shapeSelection)
{
    KIS_SAFE_ASSERT_RECOVER_RETURN_VALUE(!m_d->shapeSelection, nullptr);

    shapeSelection->setResolutionProxy(m_d->resolutionProxy);
    return new ChangeShapeSelectionCommand(this, shapeSelection);
}

KisPixelSelectionSP KisSelection::projection() const
{
    return m_d->pixelSelection;
}

void KisSelection::updateProjection(const QRect &rc)
{
    QReadLocker l(&m_d->shapeSelectionPointerLock);

    if(m_d->shapeSelection) {
        m_d->shapeSelection->renderToProjection(m_d->pixelSelection, rc);
        m_d->pixelSelection->setOutlineCache(m_d->shapeSelection->outlineCache());
    }
}

void KisSelection::updateProjection()
{
    QReadLocker l(&m_d->shapeSelectionPointerLock);

    if(m_d->shapeSelection) {
        m_d->pixelSelection->clear();
        m_d->shapeSelection->renderToProjection(m_d->pixelSelection);
        m_d->pixelSelection->setOutlineCache(m_d->shapeSelection->outlineCache());
    }
}

void KisSelection::setVisible(bool visible)
{
    bool needsNotification = visible != m_d->isVisible;

    m_d->isVisible = visible;

    if (needsNotification) {
        notifySelectionChanged();
    }
}

bool KisSelection::isVisible()
{
    return m_d->isVisible;
}

bool KisSelection::isTotallyUnselected(const QRect & r) const
{
    return m_d->pixelSelection->isTotallyUnselected(r);
}

QRect KisSelection::selectedRect() const
{
    return m_d->pixelSelection->selectedRect();
}

QRect KisSelection::selectedExactRect() const
{
    return m_d->pixelSelection->selectedExactRect();
}

qint32 KisSelection::x() const
{
    return m_d->pixelSelection->x();
}

qint32 KisSelection::y() const
{
    return m_d->pixelSelection->y();
}

void KisSelection::setX(qint32 x)
{
    QReadLocker l(&m_d->shapeSelectionPointerLock);

    Q_ASSERT(m_d->pixelSelection);

    qint32 delta = x - m_d->pixelSelection->x();
    m_d->pixelSelection->setX(x);
    if (m_d->shapeSelection) {
        m_d->shapeSelection->moveX(delta);
    }
}

void KisSelection::setY(qint32 y)
{
    QReadLocker l(&m_d->shapeSelectionPointerLock);

    Q_ASSERT(m_d->pixelSelection);

    qint32 delta = y - m_d->pixelSelection->y();
    m_d->pixelSelection->setY(y);
    if (m_d->shapeSelection) {
        m_d->shapeSelection->moveY(delta);
    }
}

void KisSelection::setDefaultBounds(KisDefaultBoundsBaseSP bounds)
{
    m_d->pixelSelection->setDefaultBounds(bounds);
}

void KisSelection::setResolutionProxy(KisImageResolutionProxySP proxy)
{
    m_d->resolutionProxy = proxy;
    if (m_d->shapeSelection) {
        m_d->shapeSelection->setResolutionProxy(proxy);
    }
}

KisImageResolutionProxySP KisSelection::resolutionProxy() const
{
    return m_d->resolutionProxy;
}

void KisSelection::clear()
{
    QReadLocker readLocker(&m_d->shapeSelectionPointerLock);

    if (m_d->shapeSelection) {
        readLocker.unlock();
        QWriteLocker writeLocker(&m_d->shapeSelectionPointerLock);
        if (m_d->shapeSelection) {
            Private::safeDeleteShapeSelection(m_d->shapeSelection, this);
            m_d->shapeSelection = 0;
        }
    }

    m_d->pixelSelection->clear();
}

KUndo2Command* KisSelection::flatten()
{
    QReadLocker readLocker(&m_d->shapeSelectionPointerLock);

    KUndo2Command *command = 0;

    if (m_d->shapeSelection) {
        command = m_d->shapeSelection->resetToEmpty();
        readLocker.unlock();

        if (command) {
            KisCommandUtils::CompositeCommand *cmd = new KisCommandUtils::CompositeCommand();
            cmd->addCommand(command);
            cmd->addCommand(new ChangeShapeSelectionCommand(this, nullptr));
            command = cmd;
        } else {
            command = new ChangeShapeSelectionCommand(this, nullptr);
        }
    }

    return command;
}

void KisSelection::notifySelectionChanged()
{
    KisNodeWSP parentNode;
    if (!(parentNode = this->parentNode())) return;

    KisNodeGraphListener *listener;
    if (!(listener = parentNode->graphListener())) return;

    listener->notifySelectionChanged();
}

void KisSelection::requestCompressedProjectionUpdate(const QRect &rc)
{
    m_d->updateCompressor->requestUpdate(rc);
}

quint8 KisSelection::selected(qint32 x, qint32 y) const
{
    KisHLineConstIteratorSP iter = m_d->pixelSelection->createHLineConstIteratorNG(x, y, 1);

    const quint8 *pix = iter->oldRawData();

    return *pix;
}

