# -*- coding: utf-8 -*-
from __future__ import (absolute_import, unicode_literals, division,
                        print_function)
                  
import functools
import contextlib
import warnings
import numpy as np
import numpy.fft

__all__ = ['complexmp', 'ignoredivide', 'remove_piston', 'circle_aperture',
    'fftgrid', 'shapegrid', 'shapestr', 'remove_tiptilt', 'apply_tiptilt']

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
    ap : array-like, (n x m)
        the grid of valid measurement positions
    phi : array-like, (n x m)
        the phase measurement
    
    Returns
    -------
    phi_np : array-like, (n x m)
        the slope measurement with piston removed.
    ps : float
        the piston value.
        
    
    Notes
    -----
    In slopes, piston is tip or tilt. Removing piston from slopes
    will have the effect of removing the tip or tilt of the slopes.
    
    """
    ps = np.sum(phi * ap) / np.sum(ap)
    phi_np = phi - (ps * ap)
    return (phi_np, ps)
    
def remove_tiptilt(ap, phi):
    """Remove tip and tilt from a phase.
    
    Parameters
    ----------
    ap : array-like, (n x m)
        the grid of valid measurement positions
    phi : array-like, (n x m)
        the phase measurement
    
    Returns
    -------
    phi_np : array-like, (n x m)
        the slope measurement with tip and tilt removed.
    tilt_x : float
        the x tilt value which was removed
    tilt_y : float
        the y tilt value which was removed
    
    """
    y, x = shapegrid(phi.shape)
    
    x = ap * ((x * ap) - (x * ap).sum() / ap.sum())
    y = ap * ((y * ap) - (y * ap).sum() / ap.sum())
    
    tip  = np.sum(phi * ap * y) / np.sum(ap * y**2)
    tilt = np.sum(phi * ap * x) / np.sum(ap * x**2)
    
    phi_ntt = phi - (tip * y * ap) - (tilt * x * ap)
    return (phi_ntt, tilt, tip)
    
def apply_tiptilt(ap, phi, tt_x, tt_y):
    """Apply tip and tilt to a phase."""
    y, x = shapegrid(phi.shape)
    
    x = ap * ((x * ap) - (x * ap).sum() / ap.sum())
    y = ap * ((y * ap) - (y * ap).sum() / ap.sum())
    
    phi_tt = phi + (tt_y * y * ap) + (tt_x * x * ap)
    return phi_tt
    
def circle_aperture(shape, r):
    """Create a circle aperture."""
    y, x = shapegrid(shape)
    return ((x**2 + y**2) < r**2).astype(np.int)
    
def fftgrid(shape, scale=1.0/(2.0*np.pi)):
    """FFT Grid"""
    ff = [np.fft.fftfreq(s, scale) for s in reversed(shape)]
    return tuple(reversed(np.meshgrid(*ff)))
    
def shapegrid(shape, centered=True):
    """A grid shaped for numpy."""
    if centered:
        ii = [ np.arange(s) - ((s-1) / 2.0) for s in reversed(shape) ]
    else:
        ii = [ np.arange(s) for s in reversed(shape) ]
    return tuple(reversed(np.meshgrid(*ii)))
    
def shapestr(shape):
    """Return a string formatted shape tuple."""
    return "({0})".format("x".join("{0:d}".format(s) for s in shape))
    
def is_hermitian(matrix):
    """docstring for is_hermitian"""
    matrix = np.asmatrix(matrix, dtype=np.complex)
    if matrix.shape != tuple(reversed(matrix.shape)):
        return False
    return np.allclose(matrix.H, matrix)
    
def make_hermitian(matrix):
    """Make an array hermitian, assuming the top-right triangle is valid."""
    matrix = np.asmatrix(matrix, dtype=np.complex)
    if matrix.shape != tuple(reversed(matrix.shape)):
        raise ValueError("Can't make a non-square matrix hermitian.")
    triu = np.matrix(np.triu(matrix))
    matrix_h = triu + triu.H
    matrix_h -= np.diag(matrix_h.diagonal()) / 2.0
    return matrix_h

def fill_hermitian_fft(array):
    """Fill an FFT-shifted array so it is hermitian."""
    array_shifted = np.fft.fftshift(np.asarray(array, dtype=np.complex))
    array_hermitian = make_hermitian(array_shifted)
    return np.fft.ifftshift(array_hermitian.view(np.ndarray))
