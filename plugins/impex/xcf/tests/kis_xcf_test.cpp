/*
 * SPDX-FileCopyrightText: 2007 Cyrille Berger <cberger@cberger.net>
 *
 *  SPDX-License-Identifier: GPL-2.0-or-later
 */

#include "kis_xcf_test.h"


#include <simpletest.h>
#include <QCoreApplication>

#include <testui.h>

#include "filestest.h"

#ifndef FILES_DATA_DIR
#error "FILES_DATA_DIR not set. A directory with the data used for testing the importing of files in krita"
#endif

const QString XcfMimetype = "image/x-xcf";


void KisXCFTest::testFiles()
{
    TestUtil::testFiles(QString(FILES_DATA_DIR) + "/sources", QStringList(), QString(), 1);
}


void KisXCFTest::testImportFromWriteonly()
{
    TestUtil::testImportFromWriteonly(XcfMimetype);
}

/*
void KisXCFTest::testExportToReadonly()
{
    TestUtil::testExportToReadonly(XcfMimetype);
}
*/


void KisXCFTest::testImportIncorrectFormat()
{
    TestUtil::testImportIncorrectFormat(XcfMimetype);
}


KISTEST_MAIN(KisXCFTest)

