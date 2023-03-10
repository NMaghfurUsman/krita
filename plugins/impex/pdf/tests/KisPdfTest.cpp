/*
 * SPDX-FileCopyrightText: 2019 Agata Cacko <cacko.azh@gmail.com>
 *
 *  SPDX-License-Identifier: GPL-2.0-or-later
 */

#include "KisPdfTest.h"


#include <simpletest.h>
#include <QCoreApplication>

#include "filestest.h"

#ifndef FILES_DATA_DIR
#error "FILES_DATA_DIR not set. A directory with the data used for testing the importing of files in krita"
#endif


const QString PdfMimetype = "image/x-gimp-brush";



void KisPdfTest::testImportFromWriteonly()
{
    TestUtil::testImportFromWriteonly(PdfMimetype);
}


void KisPdfTest::testExportToReadonly()
{
    TestUtil::testExportToReadonly(PdfMimetype);
}


void KisPdfTest::testImportIncorrectFormat()
{
    TestUtil::testImportIncorrectFormat(PdfMimetype);
}



KISTEST_MAIN(KisPdfTest)


