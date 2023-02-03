#include "NPointPerspective.h"
#include "KoPointerEvent.h"
#include "kis_painting_assistant.h"
#include <qglobal.h>

NPointPerspective::NPointPerspective(const QString &id, const QString &name)
    : KisPaintingAssistant(id, name) {}

NPointPerspective::NPointPerspective(
    const KisPaintingAssistant &rhs,
    QMap<KisPaintingAssistantHandleSP, KisPaintingAssistantHandleSP> &handleMap)
  : KisPaintingAssistant(rhs, handleMap)
{}

NPointPerspective::~NPointPerspective() {}

void NPointPerspective::realignSideHandles(KisPaintingAssistantHandleSP dragged_handle) {
    Q_UNUSED(dragged_handle);
}

void NPointPerspective::realignVanishingPoint(KisPaintingAssistantHandleSP dragged_handle, KoPointerEvent* event, QPointF* drag_start, QPointF* adjustment) {
    Q_UNUSED(dragged_handle);
    Q_UNUSED(drag_start);
    Q_UNUSED(adjustment);
}

void NPointPerspective::initSideHandles()
{

}
