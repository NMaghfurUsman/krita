/*
 *  SPDX-FileCopyrightText: 2004, 2007-2010 Cyrille Berger <cberger@cberger.net>
 *  SPDX-FileCopyrightText: 2018 Ivan Santa Maria <ghevan@gmail.com>
 *  SPDX-FileCopyrightText: 2022 L. E. Segovia <amy@amyspark.me>
 *
 *  SPDX-License-Identifier: GPL-2.0-or-later
 */

#include <cmath>


#include <QDomDocument>

#include "KoMultiArchBuildSupport.h"
#include "kis_base_mask_generator.h"
#include "kis_fast_math.h"
#include "kis_rect_mask_generator.h"
#include "kis_rect_mask_generator_p.h"


#include "kis_brush_mask_applicator_factories.h"
#include "kis_brush_mask_applicator_base.h"

#include <qnumeric.h>

KisRectangleMaskGenerator::KisRectangleMaskGenerator(qreal radius, qreal ratio, qreal fh, qreal fv, int spikes, bool antialiasEdges)
    : KisMaskGenerator(radius, ratio, fh, fv, spikes, antialiasEdges, RECTANGLE, DefaultId), d(new Private)
{
    setScale(1.0, 1.0);

    // store the variable locally to allow vector implementation read it easily
    d->copyOfAntialiasEdges = antialiasEdges;
    d->applicator.reset(createOptimizedClass<MaskApplicatorFactory<KisRectangleMaskGenerator>>(this));
}

KisRectangleMaskGenerator::KisRectangleMaskGenerator(const KisRectangleMaskGenerator &rhs)
    : KisMaskGenerator(rhs),
      d(new Private(*rhs.d))
{
    d->applicator.reset(createOptimizedClass<MaskApplicatorFactory<KisRectangleMaskGenerator>>(this));
}

KisMaskGenerator* KisRectangleMaskGenerator::clone() const
{
    return new KisRectangleMaskGenerator(*this);
}

KisRectangleMaskGenerator::~KisRectangleMaskGenerator()
{
}

void KisRectangleMaskGenerator::setScale(qreal scaleX, qreal scaleY)
{
    KisMaskGenerator::setScale(scaleX, scaleY);

    d->xcoeff = 2.0 / effectiveSrcWidth();
    d->ycoeff = 2.0 / effectiveSrcHeight();
    d->xfadecoeff = (horizontalFade() == 0) ? 1 : (2.0 / (horizontalFade() * effectiveSrcWidth()));
    d->yfadecoeff = (verticalFade() == 0)   ? 1 : (2.0 / (verticalFade() * effectiveSrcHeight()));

    setSoftness(this->softness());
}

void KisRectangleMaskGenerator::setSoftness(qreal softness)
{
    KisMaskGenerator::setSoftness(softness);
    qreal safeSoftnessCoeff = qreal(1.0) / qMax(qreal(0.01), softness);

    d->transformedFadeX = d->xfadecoeff * safeSoftnessCoeff;
    d->transformedFadeY = d->yfadecoeff * safeSoftnessCoeff;
}

bool KisRectangleMaskGenerator::shouldVectorize() const
{
    return !shouldSupersample() && spikes() == 2;
}

KisBrushMaskApplicatorBase *KisRectangleMaskGenerator::applicator() const
{
    return d->applicator.data();
}

void KisRectangleMaskGenerator::setMaskScalarApplicator()
{
    d->applicator.reset(
        createScalarClass<MaskApplicatorFactory<KisRectangleMaskGenerator>>(
            this));
}

quint8 KisRectangleMaskGenerator::valueAt(qreal x, qreal y) const
{
    if (isEmpty()) return 255;
    qreal xr = qAbs(x /*- m_xcenter*/);
    qreal yr = qAbs(y /*- m_ycenter*/);
    fixRotation(xr, yr);

    xr = qAbs(xr);
    yr = qAbs(yr);

    qreal nxr = xr * d->xcoeff;
    qreal nyr = yr * d->ycoeff;

    if (nxr > 1.0 || nyr > 1.0) return 255;

    if (antialiasEdges()) {
        xr += 1.0;
        yr += 1.0;
    }

    qreal fxr = xr * d->transformedFadeX;
    qreal fyr = yr * d->transformedFadeY;

    qreal fxnorm = nxr * (fxr - 1.0) / (fxr - nxr);
    qreal fynorm = nyr * (fyr - 1.0) / (fyr - nyr);

    qreal retValue = 0;

    if(fxr > 1.0) {
        retValue = fxnorm;
     }

    if (fxnorm < fynorm && fyr > 1.0) {
        retValue = fynorm;
     }

    return retValue * 255;
}

