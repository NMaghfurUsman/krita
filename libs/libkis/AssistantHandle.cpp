
#include <AssistantHandle.h>

struct AssistantHandle::Private {
    Private() {}
    KisPaintingAssistantHandle* handle;
};

AssistantHandle::AssistantHandle(KisPaintingAssistantHandle* handle, QObject *parent)
    : QObject(parent)
    , d(new Private)
{
  ENTER_FUNCTION();
  d->handle = handle;
}

AssistantHandle::~AssistantHandle()
{
    delete d;
}

bool AssistantHandle::operator==(const AssistantHandle &other) const
{
    return (d->handle == other.d->handle);
}

bool AssistantHandle::operator!=(const AssistantHandle &other) const
{
    return !(operator==(other));
}

qreal AssistantHandle::x()
{
  return d->handle->x();
}

qreal AssistantHandle::y()
{
  return d->handle->y();
}

void AssistantHandle::setX(qreal val)
{
  d->handle->setX(val);
  d->handle->uncache();
}

void AssistantHandle::setY(qreal val)
{
  d->handle->setY(val);
  d->handle->uncache();
}
