/* This file is part of the KDE project
 *
 * SPDX-FileCopyrightText: 2009 Inge Wallin <inge@lysator.liu.se>
 *
 * SPDX-License-Identifier: LGPL-2.0-or-later
 */

#include "ImageShapePlugin.h"

#include <kpluginfactory.h>

#include <KoShapeRegistry.h>

#include "ImageShapeFactory.h"

K_PLUGIN_FACTORY_WITH_JSON(ImageShapePluginFactory, "krita_shape_image.json", registerPlugin<ImageShapePlugin>();)

ImageShapePlugin::ImageShapePlugin(QObject *parent, const QVariantList &)
    : QObject(parent)
{
    KoShapeRegistry::instance()->add(new ImageShapeFactory());
}

#include <ImageShapePlugin.moc>
