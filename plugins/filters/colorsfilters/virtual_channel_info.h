/*
 *  SPDX-FileCopyrightText: 2015 Dmitry Kazakov <dimula73@gmail.com>
 *
 *  SPDX-License-Identifier: GPL-2.0-or-later
 */

#ifndef __VIRTUAL_CHANNEL_INFO_H
#define __VIRTUAL_CHANNEL_INFO_H

#include <KoChannelInfo.h>
class KoColorSpace;


/**
 * This class represents a virtual channel that can have a curve in
 * curves filter. Virtual channel can be of various types:
 *
 * - REAL --- represents a real color channel of the image,
 *            like R, G, B or A
 *
 * - LIGHTNESS --- lightness virtual channel: represents L channel
 *                 of the image separation into Lab.
 *
 * - ALL_COLORS --- represents a grouped channel, combining all the
 *                  color channels of the image. E.g. R+G+B of an RGB
 *                  image
 */

class VirtualChannelInfo
{
public:
    enum Type {
        REAL,
        HUE,
        SATURATION,
        LIGHTNESS,
        ALL_COLORS
    };

    VirtualChannelInfo();

    VirtualChannelInfo(Type type, int pixelIndex, KoChannelInfo *realChannelInfo, const KoColorSpace *cs);

    /**
     * \return a pointer to a KoChannelInfo structure *iff* the
     *         channel type is 'REAL'. Returns null of all the
     *         other types.
     */
    KoChannelInfo* channelInfo() const;

    /**
     * Index of this channel in a pixel.
     *
     * \return -1 for all virtual channels.
     */
    int pixelIndex() const;

    Type type() const;
    QString name() const;

    KoChannelInfo::enumChannelValueType valueType() const;
    int channelSize() const;

    bool isAlpha() const;

private:
    Type m_type {REAL};
    int m_pixelIndex {0};
    KoChannelInfo *m_realChannelInfo {nullptr};

    QString m_nameOverride;
    KoChannelInfo::enumChannelValueType m_valueTypeOverride {KoChannelInfo::FLOAT32};
    int m_channelSizeOverride {4};
};

#endif /* __VIRTUAL_CHANNEL_INFO_H */
