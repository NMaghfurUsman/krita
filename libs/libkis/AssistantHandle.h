
#ifndef LIBKIS_ASSISTANT_HANDLE_H
#define LIBKIS_ASSISTANT_HANDLE_H

#include <QObject>
#include "kritalibkis_export.h"
#include <kis_painting_assistant.h>
#include <kis_debug.h>

class AssistantHandle;

class KRITALIBKIS_EXPORT AssistantHandle : public QObject
{

  Q_OBJECT

public:
  explicit AssistantHandle(KisPaintingAssistantHandle* handle, QObject *parent = 0);
  ~AssistantHandle() override;
  
  bool operator==(const AssistantHandle &other) const;
  bool operator!=(const AssistantHandle &other) const;

public Q_SLOTS:

  qreal x();
  qreal y();
  void setX(qreal val);
  void setY(qreal val);

 private:
  struct Private;
  Private *const d;
};

#endif // LIBKIS_ASSISTANT_HANDLE_H
