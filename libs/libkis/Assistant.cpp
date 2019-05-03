
#include <Assistant.h>

struct Assistant::Private {
    Private() {}
  QSharedPointer<KisPaintingAssistant> assis;
};

Assistant::Assistant(QSharedPointer<KisPaintingAssistant> assis, QObject *parent)
    : QObject(parent)
    , d(new Private)
{
  ENTER_FUNCTION();
  d->assis = assis;
}

Assistant::~Assistant()
{
    delete d;
}

bool Assistant::operator==(const Assistant &other) const
{
    return (d->assis == other.d->assis);
}

bool Assistant::operator!=(const Assistant &other) const
{
    return !(operator==(other));
}

QString Assistant::id()
{
  if (!d->assis) return "";
  return d->assis->id();
}

QString Assistant::name()
{
  if (!d->assis) return "";
  return d->assis->name();
}

QList<AssistantHandle*> Assistant::handles() const
{
  QList<AssistantHandle*> handles = {};
  Q_FOREACH (KisSharedPtr<KisPaintingAssistantHandle> handle, d->assis->handles()){
    handles.append(new AssistantHandle(handle.data()));
  }
  return handles;
}

QSharedPointer<KisPaintingAssistant> Assistant::assistant()
{
  return d->assis;
}
