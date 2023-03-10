/*
 * SPDX-FileCopyrightText: 2019 Agata Cacko <cacko.azh@gmail.com>
 *
 *  SPDX-License-Identifier: GPL-2.0-or-later
 */

#include "KisRawTest.h"


#include <simpletest.h>
#include <QCoreApplication>

#include "filestest.h"

#ifndef FILES_DATA_DIR
#error "FILES_DATA_DIR not set. A directory with the data used for testing the importing of files in krita"
#endif


const QString RawMimetype = "image/x-krita-raw";



void KisRawTest::testImportFromWriteonly()
{
    TestUtil::testImportFromWriteonly(RawMimetype);
}


void KisRawTest::testImportIncorrectFormat()
{
    TestUtil::testImportIncorrectFormat(RawMimetype);
}



KISTEST_MAIN(KisRawTest)


