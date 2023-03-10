/*
 * SPDX-FileCopyrightText: 2019 Agata Cacko <cacko.azh@gmail.com>
 *
 *  SPDX-License-Identifier: GPL-2.0-or-later
 */

#include "KisBrushTest.h"


#include <simpletest.h>
#include <QCoreApplication>

#include "filestest.h"

#ifndef FILES_DATA_DIR
#error "FILES_DATA_DIR not set. A directory with the data used for testing the importing of files in krita"
#endif


const QString BrushMimetype = "image/x-gimp-brush";



void KisBrushTest::testImportFromWriteonly()
{
    TestUtil::testImportFromWriteonly(BrushMimetype);
}


void KisBrushTest::testExportToReadonly()
{
    TestUtil::testExportToReadonly(BrushMimetype);
}


void KisBrushTest::testImportIncorrectFormat()
{
    TestUtil::testImportIncorrectFormat(BrushMimetype);
}



KISTEST_MAIN(KisBrushTest)


