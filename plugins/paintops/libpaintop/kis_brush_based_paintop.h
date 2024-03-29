/*
 *  SPDX-FileCopyrightText: 2008 Boudewijn Rempt <boud@valdyas.org>
 *
 *  SPDX-License-Identifier: GPL-2.0-or-later
 */
#ifndef KIS_BRUSH_BASED_PAINTOP_H
#define KIS_BRUSH_BASED_PAINTOP_H

#include "kritapaintop_export.h"
#include <brushengine/kis_paintop.h>
#include "kis_dab_cache.h"
#include "kis_brush.h"
#include "kis_texture_option.h"
#include "kis_precision_option.h"
#include <kis_threaded_text_rendering_workaround.h>
#include <KisMirrorOption.h>

class KisPropertiesConfiguration;
class KisDabCache;
class KisResourcesInterface;
struct KisAirbrushOptionData;
class KisSpacingOption;

/// Internal
class TextBrushInitializationWorkaround
{
public:
    TextBrushInitializationWorkaround();
    ~TextBrushInitializationWorkaround();
    static TextBrushInitializationWorkaround* instance();

    void preinitialize(KisPaintOpSettingsSP settings);

    KisBrushSP tryGetBrush(const KisPropertiesConfigurationSP settings);


private:
    KisBrushSP m_brush;
    KisPropertiesConfigurationSP m_settings;
};


/**
 * This is a base class for paintops that use a KisBrush or derived
 * brush to paint with. This is mainly important for the spacing
 * generation.
 */
class PAINTOP_EXPORT KisBrushBasedPaintOp : public KisPaintOp
{

public:

    KisBrushBasedPaintOp(const KisPaintOpSettingsSP settings, KisPainter* painter, KisBrushTextureFlags textureFlags = None);
    ~KisBrushBasedPaintOp() override;

    bool checkSizeTooSmall(qreal scale);

    KisSpacingInformation effectiveSpacing(qreal scale) const;
    KisSpacingInformation effectiveSpacing(qreal scale, qreal rotation, const KisPaintInformation &pi) const;
    KisSpacingInformation effectiveSpacing(qreal scale, qreal rotation, const KisSpacingOption &spacingOption, const KisPaintInformation &pi) const;
    KisSpacingInformation effectiveSpacing(qreal scale, qreal rotation,
                                           const KisAirbrushOptionData *airbrushOption,
                                           const KisSpacingOption *spacingOption,
                                           const KisPaintInformation &pi) const;

    ///Reimplemented, false if brush is 0
    bool canPaint() const override;

#ifdef HAVE_THREADED_TEXT_RENDERING_WORKAROUND
    typedef int needs_preinitialization;
    static void preinitializeOpStatically(KisPaintOpSettingsSP settings);
#endif /* HAVE_THREADED_TEXT_RENDERING_WORKAROUND */

    static QList<KoResourceLoadResult> prepareLinkedResources(const KisPaintOpSettingsSP settings, KisResourcesInterfaceSP resourcesInterface);
    static QList<KoResourceLoadResult> prepareEmbeddedResources(const KisPaintOpSettingsSP settings, KisResourcesInterfaceSP resourcesInterface);

private:
    KisSpacingInformation effectiveSpacing(qreal dabWidth, qreal dabHeight, qreal extraScale, bool isotropicSpacing, qreal rotation, bool axesFlipped) const;

protected: // XXX: make private!
    KisDabCache *m_dabCache;
    KisBrushSP m_brush;

private:
    KisTextureOption m_textureOption;

protected:
    KisMirrorOption m_mirrorOption;
    KisPrecisionOption m_precisionOption;
};

#endif
