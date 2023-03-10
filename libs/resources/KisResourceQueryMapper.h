/*
 *  SPDX-License-Identifier: GPL-3.0-or-later
 */

#ifndef KISRESOURCEQUERYMAPPER_H
#define KISRESOURCEQUERYMAPPER_H

#include <QSqlQuery>
#include <QMap>

class QImage;

class KisResourceQueryMapper
{
public:
    /**
     * @brief variantFromResourceQuery returns a QVariant for the given column and or role
     * @param query the query: it supposed to be in the right position already
     */
    static QVariant variantFromResourceQuery(const QSqlQuery &query, int column, int role, bool useResourcePrefix);

private:
    static QImage getThumbnailFromQuery(const QSqlQuery &query, bool useResourcePrefix);
};

#endif // KISRESOURCEQUERYMAPPER_H
