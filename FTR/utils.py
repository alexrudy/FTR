# -*- coding: utf-8 -*-
from __future__ import (absolute_import, unicode_literals, division,
                        print_function)
                  
import functools, contextlib
import numpy as np
import numpy.fft

__all__ = ['complexmp', 'ignoredivide', 'remove_piston', 'circle_aperture',
    'fftgrid', 'shapegrid', 'shapestr']

def complexmp(mag, phase):
    """Return a complex number from the magnitude and phase of that number."""
    return mag * np.cos(phase) + mag * 1j * np.sin(phase)
    
@contextlib.contextmanager
def ignoredivide():
    """A context manager for ignoring division"""
    errsettings = np.seterr(all='ignore')
    yield
    np.seterr(**errsettings)
    
def remove_piston(ap, phi):
    """Remove piston.
    
    Piston is the average position of the array in valid apertures.
    
    Parameters
    ----------
    ap, array-like, (n x m):
        the grid of valid measurement positions
    phi, array-like, (n x m):
        the phase measurement
    
    Returns
    -------
    phi_np, array-like, (n x m):
        the slope measurement with piston removed.
    ps, float:
        the tip/tilt value.
        
    
    Notes
    -----
    In slopes, piston is tip or tilt. Removing piston from slopes
    will have the effect of removing the tip or tilt of the slopes.
    
    """
    ps = np.sum(phi * ap) / np.sum(ap)
    phi_np = phi - (ps * ap)
    return (phi_np, ps)
    
def circle_aperture(shape, r):
    """Create a circle aperture."""
    y, x = shapegrid(shape)
    return ((x**2 + y**2) <= r**2).astype(np.int)
    
def fftgrid(shape, scale=1.0/(2.0*np.pi)):
    """FFT Grid"""
    ff = [np.fft.fftfreq(s, scale) for s in reversed(shape)]
    return tuple(reversed(np.meshgrid(*ff)))
    
def shapegrid(shape, centered=True):
    """A grid shaped for numpy."""
    if centered:
        ii = [ np.arange(s) - (s / 2.0) for s in reversed(shape) ]
    else:
        ii = [ np.arange(s) for s in reversed(shape) ]
    return tuple(reversed(np.meshgrid(*ii)))
    
def shapestr(shape):
    """Return a string formatted shape tuple."""
    return "({0})".format("x".join("{0:d}".format(s) for s in shape))