/* This file is part of the KDE project
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Library General Public
 * License as published by the Free Software Foundation; either
 * version 2 of the License, or (at your option) any later version.
 *
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Library General Public License for more details.
 *
 * You should have received a copy of the GNU Library General Public License
 * along with this library; see the file COPYING.LIB.  If not, write to
 * the Free Software Foundation, Inc., 51 Franklin Street, Fifth Floor,
 * Boston, MA 02110-1301, USA.
 */

#ifndef TREESHAPE_H
#define TREESHAPE_H

#include "KoShapeContainer.h"


#define TREESHAPEID "TreeShape"

class KoViewConverter;
class KoShapeSavingContext;
class KoShapeLoadingContext;
class KoShape;
class KoConnectionShape;
class Layout;
class QPainter;

class TreeShape : public KoShapeContainer
{

public:
    TreeShape();
    TreeShape(KoShape *root);
    virtual ~TreeShape();
    virtual void addChild(KoShape *tree, KoShape *connector);
    virtual KoShape* connector(KoShape *shape);
    virtual QList<KoShape*> addNewChild();
    virtual void setNextShape(KoShape *shape);
    virtual KoShape* nextShape();
    virtual void setZIndex(int zIndex);
    virtual void paintComponent(QPainter &painter, const KoViewConverter &converter);
    virtual bool hitTest(const QPointF &position) const;
    virtual void saveOdf(KoShapeSavingContext &context) const;
    virtual bool loadOdf(const KoXmlElement &element, KoShapeLoadingContext &context);

private:
//     virtual void shapeChanged(ChangeType type, KoShape *shape = 0);
    virtual Layout *layout() const;

    KoShape *m_nextShape;
    uint m_structure;
};

#endif // KOTREESHAPE_H
