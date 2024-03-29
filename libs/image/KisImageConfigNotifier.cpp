/*
 *  SPDX-FileCopyrightText: 2017 Dmitry Kazakov <dimula73@gmail.com>
 *
 *  SPDX-License-Identifier: GPL-2.0-or-later
 */

#include "KisImageConfigNotifier.h"

#include <QGlobalStatic>

#include <kis_debug.h>
#include "kis_signal_compressor.h"

Q_GLOBAL_STATIC(KisImageConfigNotifier, s_instance)

struct KisImageConfigNotifier::Private
{
    Private()
        : updateCompressor(300, KisSignalCompressor::FIRST_ACTIVE)
        , autoKeyframeUpdateCompressor(300, KisSignalCompressor::FIRST_ACTIVE)
    {}

    KisSignalCompressor updateCompressor;
    KisSignalCompressor autoKeyframeUpdateCompressor;
};

KisImageConfigNotifier::KisImageConfigNotifier()
    : m_d(new Private)
{
    connect(&m_d->updateCompressor, SIGNAL(timeout()), SIGNAL(configChanged()));
    connect(&m_d->updateCompressor, SIGNAL(timeout()), SIGNAL(autoKeyFrameConfigurationChanged()));
    connect(&m_d->autoKeyframeUpdateCompressor, SIGNAL(timeout()), SIGNAL(autoKeyFrameConfigurationChanged()));
}

KisImageConfigNotifier::~KisImageConfigNotifier()
{
}

KisImageConfigNotifier *KisImageConfigNotifier::instance()
{
    return s_instance;
}

void KisImageConfigNotifier::notifyConfigChanged()
{
    m_d->updateCompressor.start();
}

void KisImageConfigNotifier::notifyAutoKeyFrameConfigurationChanged()
{
    m_d->autoKeyframeUpdateCompressor.start();
}
