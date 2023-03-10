/*
 * SPDX-FileCopyrightText: 2019 Agata Cacko <cacko.azh@gmail.com>
 *
 *  SPDX-License-Identifier: GPL-2.0-or-later
 */

#include "KisTgaTest.h"


#include <simpletest.h>
#include <QCoreApplication>

#include "filestest.h"

#ifndef FILES_DATA_DIR
#error "FILES_DATA_DIR not set. A directory with the data used for testing the importing of files in krita"
#endif


const QString TgaMimetype = "image/x-gimp-brush";



void KisTgaTest::testImportFromWriteonly()
{
    TestUtil::testImportFromWriteonly(TgaMimetype);
}


void KisTgaTest::testExportToReadonly()
{
    TestUtil::testExportToReadonly(TgaMimetype);
}


void KisTgaTest::testImportIncorrectFormat()
{
    TestUtil::testImportIncorrectFormat(TgaMimetype);
}



KISTEST_MAIN(KisTgaTest)


