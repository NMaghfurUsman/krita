/*
 * This file is part of the KDE project
 * SPDX-FileCopyrightText: 2013 Arjen Hiemstra <ahiemstra@heimr.nl>
 *
 * SPDX-License-Identifier: GPL-2.0-or-later
 */

#ifndef KISINPUTCONFIGURATIONPAGEITEM_H
#define KISINPUTCONFIGURATIONPAGEITEM_H

#include <QWidget>

namespace Ui
{
class KisInputConfigurationPageItem;
}

class KisActionShortcutsModel;
class KisAbstractInputAction;
/**
 * \brief A collapsible widget displaying an action, its description and associated shortcuts.
 *
 * This is used in KisInputConfigurationPage to display a list of actions and the associated
 * shortcuts, depending on the current profile.
 */
class KisInputConfigurationPageItem : public QWidget
{
    Q_OBJECT
public:
    KisInputConfigurationPageItem(QWidget *parent = nullptr, Qt::WindowFlags f = Qt::WindowFlags());
    ~KisInputConfigurationPageItem() override;

    void setAction(KisAbstractInputAction *action);

    void setWarningEnabled(bool enabled, QString additionalToolTipText = QString());

Q_SIGNALS:
    void inputConfigurationChanged();

public Q_SLOTS:
    void setExpanded(bool expand);

private Q_SLOTS:
    void deleteShortcut();

private:
    Ui::KisInputConfigurationPageItem *ui {nullptr};
    KisAbstractInputAction *m_action {nullptr};
    KisActionShortcutsModel *m_shortcutsModel {nullptr};

    const QString m_defaultToolTipText;
};

#endif // KISINPUTCONFIGURATIONPAGEITEM_H
