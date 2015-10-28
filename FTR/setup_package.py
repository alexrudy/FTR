# -*- coding: utf-8 -*-
from __future__ import absolute_import

import os
from setup_helpers import get_external_libraries, get_libFTR_extensions, get_libFTR_include_directory

def get_extensions():
    """Get the libFTR extensions"""
    include_dir = get_libFTR_include_directory()
    return get_libFTR_extensions(__file__, __name__, "*.pyx", sources=[os.path.join(include_dir,"ftr.c")])