/*
 *  SPDX-FileCopyrightText: 2019 Anna Medonosova <anna.medonosova@gmail.com>
 *
 *  SPDX-License-Identifier: GPL-2.0-or-later
 */


#ifndef KISMIRRORAXISCONFIG_H
#define KISMIRRORAXISCONFIG_H

#include <QScopedPointer>

#include "kritaui_export.h"
#include <boost/operators.hpp>

class QDomElement;
class QDomDocument;

/**
 * @brief The KisMirrorAxisConfig class stores configuration for the KisMirrorAxis
 * canvas decoration. Contents are saved to/loaded from KRA documents.
 */

class KRITAUI_EXPORT KisMirrorAxisConfig : public QObject, boost::equality_comparable<KisMirrorAxisConfig>
{
    Q_OBJECT

public:
    KisMirrorAxisConfig();
    ~KisMirrorAxisConfig();

    KisMirrorAxisConfig(const KisMirrorAxisConfig &rhs);
    KisMirrorAxisConfig& operator=(const KisMirrorAxisConfig& rhs);
    bool operator==(const KisMirrorAxisConfig& rhs) const;

    bool mirrorHorizontal() const;
    void setMirrorHorizontal(bool state);

    bool mirrorVertical() const;
    void setMirrorVertical(bool state);

    bool lockHorizontal() const;
    void setLockHorizontal(bool state);

    bool lockVertical() const;
    void setLockVertical(bool state);

    bool hideVerticalDecoration() const;
    void setHideVerticalDecoration(bool state);

    bool hideHorizontalDecoration() const;
    void setHideHorizontalDecoration(bool state);

    float handleSize() const;
    void setHandleSize(float size);

    float horizontalHandlePosition() const;
    void setHorizontalHandlePosition(float position);

    float verticalHandlePosition() const;
    void setVerticalHandlePosition(float position);

    QPointF axisPosition() const;
    void setAxisPosition(QPointF position);

    /**
     * @brief saveToXml() function for KisKraSaver
     * @param doc
     * @param tag
     * @return
     */
    QDomElement saveToXml(QDomDocument& doc, const QString &tag) const;

    /**
     * @brief loadFromXml() function for KisKraLoader
     * @param parent element
     * @return
     */
    bool loadFromXml(const QDomElement &parent);

    /**
     * @brief Check whether the config object was changed, or is the class default.
     * @return true, if the object is default; false, if the config was changed
     */
    bool isDefault() const;

private:
    class Private;
    const QScopedPointer<Private> d;
};

#endif // KISMIRRORAXISCONFIG_H
