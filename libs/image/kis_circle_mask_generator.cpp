/*
 *  SPDX-FileCopyrightText: 2004, 2007-2009 Cyrille Berger <cberger@cberger.net>
 *  SPDX-FileCopyrightText: 2010 Lukáš Tvrdý <lukast.dev@gmail.com>
 *  SPDX-FileCopyrightText: 2012 Sven Langkamp <sven.langkamp@gmail.com>
 *
 *  SPDX-License-Identifier: GPL-2.0-or-later
 */

#include <cmath>

#include <QDomDocument>

#include "kis_fast_math.h"
#include "kis_circle_mask_generator.h"
#include "kis_circle_mask_generator_p.h"
#include "kis_base_mask_generator.h"
#include "kis_brush_mask_applicator_factories.h"
#include "kis_brush_mask_applicator_base.h"


KisCircleMaskGenerator::KisCircleMaskGenerator(qreal diameter, qreal ratio, qreal fh, qreal fv, int spikes, bool antialiasEdges)
    : KisMaskGenerator(diameter, ratio, fh, fv, spikes, antialiasEdges, CIRCLE, DefaultId),
      d(new Private)
{
    setScale(1.0, 1.0);

    // store the variable locally to allow vector implementation read it easily
    d->copyOfAntialiasEdges = antialiasEdges;

    d->applicator.reset(createOptimizedClass<MaskApplicatorFactory<KisCircleMaskGenerator> >(this));
}

KisCircleMaskGenerator::KisCircleMaskGenerator(const KisCircleMaskGenerator &rhs)
    : KisMaskGenerator(rhs),
      d(new Private(*rhs.d))
{
    d->applicator.reset(createOptimizedClass<MaskApplicatorFactory<KisCircleMaskGenerator> >(this));
}

KisMaskGenerator* KisCircleMaskGenerator::clone() const
{
    return new KisCircleMaskGenerator(*this);
}

void KisCircleMaskGenerator::setScale(qreal scaleX, qreal scaleY)
{
    KisMaskGenerator::setScale(scaleX, scaleY);

    d->xcoef = 2.0 / effectiveSrcWidth();
    d->ycoef = 2.0 / effectiveSrcHeight();
    d->xfadecoef = qFuzzyCompare(horizontalFade(), 0) ? 1 : (2.0 / (horizontalFade() * effectiveSrcWidth()));
    d->yfadecoef = qFuzzyCompare(verticalFade()  , 0) ? 1 : (2.0 / (verticalFade() * effectiveSrcHeight()));
    d->transformedFadeX = d->xfadecoef * d->safeSoftnessCoeff;
    d->transformedFadeY = d->yfadecoef * d->safeSoftnessCoeff;
}

KisCircleMaskGenerator::~KisCircleMaskGenerator()
{
}

bool KisCircleMaskGenerator::shouldVectorize() const
{
    return !shouldSupersample() && spikes() == 2;
}

KisBrushMaskApplicatorBase *KisCircleMaskGenerator::applicator() const
{
    return d->applicator.data();
}

quint8 KisCircleMaskGenerator::valueAt(qreal x, qreal y) const
{
    if (isEmpty()) return 255;
    qreal xr = (x /*- m_xcenter*/);
    qreal yr = qAbs(y /*- m_ycenter*/);
    fixRotation(xr, yr);

    qreal n = norme(xr * d->xcoef, yr * d->ycoef);
    if (n > 1.0) return 255;

    // we add +1.0 to ensure correct antialiasing on the border
    if (antialiasEdges()) {
        xr = qAbs(xr) + 1.0;
        yr = qAbs(yr) + 1.0;
    }

    qreal nf = norme(xr * d->transformedFadeX,
                     yr * d->transformedFadeY);

    if (nf < 1.0) return 0;
    return 255 * n * (nf - 1.0) / (nf - n);
}

void KisCircleMaskGenerator::setSoftness(qreal softness)
{
    KisMaskGenerator::setSoftness(softness);
    d->safeSoftnessCoeff = qreal(1.0) / qMax(qreal(0.01), softness);

    d->transformedFadeX = d->xfadecoef * d->safeSoftnessCoeff;
    d->transformedFadeY = d->yfadecoef * d->safeSoftnessCoeff;
}

void KisCircleMaskGenerator::setMaskScalarApplicator()
{
    d->applicator.reset(
        createScalarClass<MaskApplicatorFactory<KisCircleMaskGenerator>>(this));
}
