
#ifndef LIBKIS_ASSISTANT_H
#define LIBKIS_ASSISTANT_H

#include <QObject>
#include "kritalibkis_export.h"
#include <kis_painting_assistant.h>
#include <kis_debug.h>
#include <AssistantHandle.h>
#include <kis_pointer_utils.h>
#include "libkis.h"

class Assistant;


class KRITALIBKIS_EXPORT Assistant : public QObject
{

  Q_OBJECT

public:
  explicit Assistant(QSharedPointer<KisPaintingAssistant> assis, QObject *parent = 0);
  ~Assistant() override;
  
  bool operator==(const Assistant &other) const;
  bool operator!=(const Assistant &other) const;

public Q_SLOTS:

  QString id();
  QString name();
  QList<AssistantHandle*> handles() const;

  QSharedPointer<KisPaintingAssistant> assistant(); // only for C++ side of API

 private:
  struct Private;
  Private *const d;
  friend class Canvas;
  friend class VPAssistant;
};

#endif // LIBKIS_ASSISTANT_H
