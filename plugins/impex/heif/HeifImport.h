/*
 *  SPDX-FileCopyrightText: 2018 Dirk Farin <farin@struktur.de>
 *  SPDX-FileCopyrightText: 2020-2021 Wolthera van Hövell tot Westerflier <griffinvalley@gmail.com>
 *  SPDX-FileCopyrightText: 2021 Daniel Novomesky <dnovomesky@gmail.com>
 *  SPDX-FileCopyrightText: 2021 L. E. Segovia <amy@amyspark.me>
 *
 *  SPDX-License-Identifier: GPL-2.0-or-later
 */

#ifndef HEIF_IMPORT_H_
#define HEIF_IMPORT_H_

#include <QVariant>

#include <KisImportExportFilter.h>

class HeifImport : public KisImportExportFilter
{
    Q_OBJECT
public:
    HeifImport(QObject *parent, const QVariantList &);
    ~HeifImport() override;
    bool supportsIO() const override { return true; }

    KisImportExportErrorCode
    convert(KisDocument *document,
            QIODevice *io,
            KisPropertiesConfigurationSP configuration = 0) override;
};

#endif
