// This file is part of PyKrita, Krita' Python scripting plugin.
//
// SPDX-FileCopyrightText: 2006 Paul Giannaros <paul@giannaros.org>
// SPDX-FileCopyrightText: 2012, 2013 Shaheed Haque <srhaque@theiet.org>
// SPDX-FileCopyrightText: 2013 Alex Turbov <i.zaufi@gmail.com>
// SPDX-FileCopyrightText: 2021 L. E. Segovia <amy@amyspark.me>
//
// SPDX-License-Identifier: LGPL-2.1-only OR LGPL-3.0-only OR LicenseRef-KDE-Accepted-LGPL
//

#include "PykritaModule.h"

#include "kis_debug.h"

struct module_state {
    PyObject *error;
};

#define GETSTATE(m) ((struct module_state*)PyModule_GetState(m))

/// \note Namespace name written in uppercase intentionally!
/// It will appear in debug output from Python plugins...
namespace PYKRITA
{
    PyObject* debug(PyObject* /*self*/, PyObject* args)
    {
        const char* text;

        if (PyArg_ParseTuple(args, "s", &text))
            dbgScript << text;
        Py_INCREF(Py_None);
        return Py_None;
    }
}                                                           // namespace PYKRITA

namespace
{
    PyMethodDef pykritaMethods[] = {
        {
            "qDebug"
            , &PYKRITA::debug
            , METH_VARARGS
            , "True KDE way to show debug info"
        }
        , { 0, 0, 0, 0 }
    };
}                                                           // anonymous namespace

//BEGIN Python module registration
static struct PyModuleDef moduledef = {
    PyModuleDef_HEAD_INIT
    , "pykrita"
    , "The pykrita module"
    , -1
    , pykritaMethods
    , 0
    , 0
    , 0
    , 0
};

#define INITERROR return NULL

PyMODINIT_FUNC PYKRITA_INIT()
{
    PyObject *pykritaModule = PyModule_Create(&moduledef);

    if (pykritaModule == NULL)
        INITERROR;

    PyModule_AddStringConstant(pykritaModule, "__file__", __FILE__);

    return pykritaModule;
}

//END Python module registration

// krita: space-indent on; indent-width 4;
#undef PYKRITA_INIT
