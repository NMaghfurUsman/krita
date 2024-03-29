/*
 *  SPDX-FileCopyrightText: 2002 Patrick Julien <freak@codepimps.org>
 *  SPDX-FileCopyrightText: 2004-2008 Boudewijn Rempt <boud@valdyas.org>
 *  SPDX-FileCopyrightText: 2004 Clarence Dang <dang@kde.org>
 *  SPDX-FileCopyrightText: 2004 Adrian Page <adrian@pagenet.plus.com>
 *  SPDX-FileCopyrightText: 2004 Cyrille Berger <cberger@cberger.net>
 *
 *  SPDX-License-Identifier: GPL-2.0-or-later
 */

#include "kis_filterop_settings.h"

#include <QDomDocument>

#include <KisFilterOptionData.h>
#include <filter/kis_filter.h>
#include <filter/kis_filter_registry.h>
#include <filter/kis_filter_configuration.h>
#include <kis_node.h>
#include <kis_image.h>
#include <kis_types.h>
#include <kis_paint_device.h>

KisFilterOpSettings::KisFilterOpSettings(KisResourcesInterfaceSP resourcesInterface)
    : KisBrushBasedPaintOpSettings(resourcesInterface)
{
    setPropertyNotSaved(KisFilterOptionData::filterConfigTag());
}

KisFilterOpSettings::~KisFilterOpSettings()
{
}

bool KisFilterOpSettings::paintIncremental()
{
    return true; // We always paint on the existing data
}

KisFilterConfigurationSP KisFilterOpSettings::filterConfig() const
{
    if (hasProperty(KisFilterOptionData::filterIdTag())) {
        KisFilterSP filter =
            KisFilterRegistry::instance()->get(
                getString(KisFilterOptionData::filterIdTag()));
        if (filter) {
            KisFilterConfigurationSP configuration = filter->factoryConfiguration(resourcesInterface());
            configuration->fromXML(getString(KisFilterOptionData::filterConfigTag()));
            return configuration;
        }
    }
    return 0;
}

void KisFilterOpSettings::toXML(QDomDocument& doc, QDomElement& root) const
{
    KisPaintOpSettings::toXML(doc, root);

    KisFilterConfigurationSP configuration = filterConfig();
    if (configuration) {
        QDomElement e = doc.createElement("filterconfig");
        configuration->toXML(doc, e);
        root.appendChild(e);
    }
}

void KisFilterOpSettings::fromXML(const QDomElement& e)
{
    KisPaintOpSettings::fromXML(e);
    QDomElement element = e.firstChildElement("filterconfig");

    if (hasProperty(KisFilterOptionData::filterIdTag())) {
        KisFilterSP filter =
            KisFilterRegistry::instance()->get(
                getString(KisFilterOptionData::filterIdTag()));
        if (filter) {
            KisFilterConfigurationSP configuration = filter->factoryConfiguration(resourcesInterface());
            configuration->fromXML(element);
            setProperty(KisFilterOptionData::filterConfigTag(), configuration->toXML());
        }
    }
}

bool KisFilterOpSettings::hasPatternSettings() const
{
    return false;
}
