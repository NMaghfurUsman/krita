/*
 *  SPDX-FileCopyrightText: 2015 Dmitry Kazakov <dimula73@gmail.com>
 *
 *  SPDX-License-Identifier: GPL-2.0-or-later
 */

#ifndef __KIS_ASL_CALLBACK_OBJECT_CATCHER_H
#define __KIS_ASL_CALLBACK_OBJECT_CATCHER_H

#include "kis_asl_object_catcher.h"

#include <QScopedPointer>
#include <functional>

#include <resources/KoAbstractGradient.h>

#include "kritapsdutils_export.h"

class KoPattern;

using ASLCallbackDouble = std::function<void(double)>;
using ASLCallbackInteger = std::function<void(int)>;
using ASLCallbackString = std::function<void(const QString &)>;
using ASLCallbackBoolean = std::function<void(bool)>;
using ASLCallbackColor = std::function<void(const KoColor &)>;
using ASLCallbackPoint = std::function<void(const QPointF &)>;
using ASLCallbackCurve = std::function<void(const QString &, const QVector<QPointF> &)>;
using ASLCallbackPattern = std::function<void(const KoPatternSP, const QString &)>;
using ASLCallbackPatternRef = std::function<void(const QString &, const QString &)>;
using ASLCallbackGradient = std::function<void(KoAbstractGradientSP)>;
using ASLCallbackNewStyle = std::function<void()>;

class KRITAPSDUTILS_EXPORT KisAslCallbackObjectCatcher : public KisAslObjectCatcher
{
public:
    KisAslCallbackObjectCatcher();
    ~KisAslCallbackObjectCatcher() override;

    void addDouble(const QString &path, double value) override;
    void addInteger(const QString &path, int value) override;
    void addEnum(const QString &path, const QString &typeId, const QString &value) override;
    void addUnitFloat(const QString &path, const QString &unit, double value) override;
    void addText(const QString &path, const QString &value) override;
    void addBoolean(const QString &path, bool value) override;
    void addColor(const QString &path, const KoColor &value) override;
    void addPoint(const QString &path, const QPointF &value) override;
    void addCurve(const QString &path, const QString &name, const QVector<QPointF> &points) override;
    void addPattern(const QString &path, const KoPatternSP pattern, const QString &patternUuid) override;
    void addPatternRef(const QString &path, const QString &patternUuid, const QString &patternName) override;
    void addGradient(const QString &path, KoAbstractGradientSP gradient) override;
    void newStyleStarted() override;

    void subscribeDouble(const QString &path, ASLCallbackDouble callback);
    void subscribeInteger(const QString &path, ASLCallbackInteger callback);
    void subscribeEnum(const QString &path, const QString &typeId, ASLCallbackString callback);
    void subscribeUnitFloat(const QString &path, const QString &unit, ASLCallbackDouble callback);
    void subscribeText(const QString &path, ASLCallbackString callback);
    void subscribeBoolean(const QString &path, ASLCallbackBoolean callback);
    void subscribeColor(const QString &path, ASLCallbackColor callback);
    void subscribePoint(const QString &path, ASLCallbackPoint callback);
    void subscribeCurve(const QString &path, ASLCallbackCurve callback);
    void subscribePattern(const QString &path, ASLCallbackPattern callback);
    void subscribePatternRef(const QString &path, ASLCallbackPatternRef callback);
    void subscribeGradient(const QString &path, ASLCallbackGradient callback);
    void subscribeNewStyleStarted(ASLCallbackNewStyle callback);

private:
    struct Private;
    const QScopedPointer<Private> m_d;
};

#endif /* __KIS_ASL_CALLBACK_OBJECT_CATCHER_H */
