# -*- coding: utf-8 -*-
from __future__ import (absolute_import, unicode_literals, division,
                        print_function)
                  
import functools, contextlib
import numpy as np


def complexmp(mag, phase):
    """Return a complex number from the magnitude and phase of that number."""
    return mag * np.cos(phase) + mag * 1j * np.sin(phase)
    
@contextlib.contextmanager
def ignoredivide():
    """A context manager for ignoring division"""
    errsettings = np.seterr(all='ignore')
    yield
    np.seterr(**errsettings)
    
def remove_tilt(ap, sl):
    """Remove tip or tilt from slopes."""
    tt = np.sum(sl * ap) / np.sum(ap)
    sl_nt = sl - tt * ap
    return (sl_nt, tt)