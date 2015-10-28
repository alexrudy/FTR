# -*- coding: utf-8 -*-
from __future__ import absolute_import

from setup_helpers import get_external_libraries, get_libFTR_extensions

def get_extensions():
    """Get the libFTR extensions"""
    return get_libFTR_extensions(__file__, __name__, "*.pyx")