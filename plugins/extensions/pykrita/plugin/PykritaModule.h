/*
 * This file is part of PyKrita, Krita' Python scripting plugin.
 *
 * SPDX-FileCopyrightText: 2013 Alex Turbov <i.zaufi@gmail.com>
 * SPDX-FileCopyrightText: 2014-2016 Boudewijn Rempt <boud@valdyas.org>
 * SPDX-FileCopyrightText: 2021 L. E. Segovia <amy@amyspark.me>
 *
 * SPDX-License-Identifier: LGPL-2.0-only OR LGPL-3.0-only
 */

#ifndef __PYKRITA_MODULE_H__
#define  __PYKRITA_MODULE_H__

#include <Python.h>
#include "config.h"

#if SIP_VERSION >= 0x0500000
#define PYKRITA_INIT PyInit_krita
#else
#define PYKRITA_INIT PyInit_pykrita
#endif

/**
 * Initializer for the built-in Python module.
 */
PyMODINIT_FUNC PYKRITA_INIT();

#endif
