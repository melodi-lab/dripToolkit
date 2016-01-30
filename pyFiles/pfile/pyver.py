#!/usr/bin/env python

import sys
import platform

if __name__ == "__main__":
    version = platform.python_version_tuple()
    if not version[0] == '2':
        raise Exception('You really need to do this under Python 2.x')

    verstring = '%s.%s' % (version[0], version[1])
    print verstring

    #'/usr/nikola/pkgs/python/include/python%s' % verstring
