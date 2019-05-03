class AssistantHandle : QObject
{
%TypeHeaderCode
#include "AssistantHandle.h"
%End
AssistantHandle(const AssistantHandle & __0);
public:
    bool operator==(const AssistantHandle &other) const;
    bool operator!=(const AssistantHandle &other) const;
public Q_SLOTS:
    qreal x();
    qreal y();
    void setX(qreal val);
    void setY(qreal val);
private:
};