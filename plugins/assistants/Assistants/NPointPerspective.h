#include "KoPointerEvent.h"
#include "kis_painting_assistant.h"

#ifndef _N_POINT_ASSISTANT_H_
#define _N_POINT_ASSISTANT_H_

class NPointPerspective : public KisPaintingAssistant {
  public:
    NPointPerspective(const QString &id, const QString &name);
    NPointPerspective(const KisPaintingAssistant &rhs,
                      QMap<KisPaintingAssistantHandleSP, KisPaintingAssistantHandleSP> &handleMap);
    virtual ~NPointPerspective();
    void virtual realignSideHandles(KisPaintingAssistantHandleSP dragged_handle);
    void virtual realignVanishingPoint(KisPaintingAssistantHandleSP dragged_handle, KoPointerEvent* event, QPointF* drag_start, QPointF* adjustment);
    void virtual initSideHandles();
};

#endif
