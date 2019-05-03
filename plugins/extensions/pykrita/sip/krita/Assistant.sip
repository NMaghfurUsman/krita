class Assistant : QObject /Abstract/
{
%TypeHeaderCode
#include "Assistant.h"
%End
Assistant(const Assistant & __0);
public:
    bool operator==(const Assistant &other) const;
    bool operator!=(const Assistant &other) const;
    virtual ~Assistant();
public Q_SLOTS:
    QString id();
    QString name();
    QList<AssistantHandle*> handles() const /Factory/;
private:
};