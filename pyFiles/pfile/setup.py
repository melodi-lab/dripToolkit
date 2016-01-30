#!/usr/bin/env python
#
# Written by John Halloran <halloj3@ee.washington.edu>
#
# Copyright (C) 2016 John Halloran
# Licensed under the Open Software License version 3.0
# See COPYING or http://opensource.org/licenses/OSL-3.0

"""
setup.py file for SWIG pFile wrapper
"""

from distutils.core import setup, Extension

setup(
    ext_modules = [
        Extension("_libpfile", sources=['error.cc', 'general.cc', 'rand.cc', 'pfile.cc', 'vbyteswapping.cc', 'libpfile_wrap.cxx'])
    ]
)
