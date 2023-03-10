/*
 *  SPDX-FileCopyrightText: 2016 Boudewijn Rempt <boud@valdyas.org>
 *
 *  SPDX-License-Identifier: LGPL-2.0-or-later
 */
#ifndef LIBKIS_RESOURCE_H
#define LIBKIS_RESOURCE_H

#include <QObject>
#include <QScopedPointer>
#include <kis_types.h>
#include "kritalibkis_export.h"
#include "libkis.h"
#include <KoResource.h>

/**
 * A Resource represents a gradient, pattern, brush tip, brush preset, palette or 
 * workspace definition.
 * 
 * @code
 * allPresets = Application.resources("preset")
 * for preset in allPresets:
 *     print(preset.name())
 * @endcode
 * 
 * Resources are identified by their type, name and filename. If you want to change
 * the contents of a resource, you should read its data using data(), parse it and
 * write the changed contents back.
 */
class KRITALIBKIS_EXPORT Resource : public QObject
{
    Q_OBJECT

public:
    Resource(int resourceId, const QString &type, const QString &name, const QString &filename, const QImage &image, QObject *parent = 0);
    Resource(KoResourceSP resource, const QString &type, QObject *parent = 0);
    ~Resource() override;
    Resource(const Resource &rhs);

    bool operator==(const Resource &other) const;
    bool operator!=(const Resource &other) const;
    Resource operator=(const Resource &rhs);


public Q_SLOTS:
    
    /**
     * Return the type of this resource. Valid types are:
     * <ul>
     * <li>pattern: a raster image representing a pattern
     * <li>gradient: a gradient
     * <li>brush: a brush tip
     * <li>preset: a brush preset
     * <li>palette: a color set
     * <li>workspace: a workspace definition.
     * </ul>
     */
    QString type() const;

    /**
     * The user-visible name of the resource.
     */
    QString name() const;
    
    /**
     * setName changes the user-visible name of the current resource.
     */
    void setName(QString value);

    /**
     * The filename of the resource, if present. Not all resources
     * are loaded from files.
     */
    QString filename() const;

    /**
     * An image that can be used to represent the resource in the
     * user interface. For some resources, like patterns, the 
     * image is identical to the resource, for others it's a mere
     * icon.
     */
    QImage image() const;
    
    /**
     * Change the image for this resource.
     */
    void setImage(QImage image);

private:

    friend class PresetChooser;
    friend class View;
    friend class Palette;
    friend class Preset;
    KoResourceSP resource() const;

    struct Private;
    QScopedPointer<Private> d;

};

#endif // LIBKIS_RESOURCE_H
