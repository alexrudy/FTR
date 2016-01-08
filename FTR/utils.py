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
    """Apply tip and tilt to a phase.
    
    Parameters
    ----------
    ap : array_like
        Aperture of valid phase points.
    phi : array_like
        Phase grid.
    tt_x : float
        X tip
    tt_y : float
        Y tip
    
    Returns
    -------
    phi_tt : array_like
        The phase, with tip and tilt terms added.
    
    
    """
    y, x = shapegrid(phi.shape)
    
    x = ap * ((x * ap) - (x * ap).sum() / ap.sum())
    y = ap * ((y * ap) - (y * ap).sum() / ap.sum())
    
    phi_tt = phi + (tt_y * y * ap) + (tt_x * x * ap)
    return phi_tt
    
def circle_aperture(shape, r):
    """Create a circle aperture.
    
    The aperture is a integer 2-d numpy array,
    with ones where the aperture is valid.
    
    Parameters
    ----------
    shape : tuple of ints
        The shape of the full aperture.
    r : float
        Radius of the inner circle within the aperture.
    
    Returns
    -------
    ap : array_like
        An integer numpy array, 1 where there is a valid
        measurement point, and 0 where it is invalid.
    
    """
    y, x = shapegrid(shape)
    return ((x**2 + y**2) < r**2).astype(np.int)
    
def fftgrid(shape, scale=1.0/(2.0*np.pi)):
    """Make a Fourier grid of x and y frequencies.
    
    This makes a grid of f_x and f_y points, by default in
    angular frequency.
    
    Parameters
    ----------
    shape : tuple of ints
        The desired shape of the Fourier grid.
    scale : float, optional
        The sampling scale for the FFT frequencies. Defaults to :math:`1/2\pi`,
        the default Fourier frequency sampling in angular frequency units.
        
    Returns
    -------
    fy : array_like
        The y Fourier frequencies.
    fx : array_like
        The x Fourier frequencies.
    
    """
    ff = [np.fft.fftfreq(s, scale) for s in reversed(shape)]
    return tuple(reversed(np.meshgrid(*ff)))
    
def shapegrid(shape, centered=True):
    """Make a centered index grid of a specific shape.
    
    Parameters
    ----------
    shape : tuple of ints
        The desired shape of the grid.
    centered : boolean, optional
        Determine whether the grid is centered.
    
    Returns
    -------
    yy : array_like
        The y indicies
    xx : array_like
        The x indicies
    
    """
    if centered:
        ii = [ np.arange(s) - ((s-1) / 2.0) for s in reversed(shape) ]
    else:
        ii = [ np.arange(s) for s in reversed(shape) ]
    return tuple(reversed(np.meshgrid(*ii)))
    
def shapestr(shape):
    """Return a string formatted shape tuple."""
    return "({0})".format("x".join("{0:d}".format(s) for s in shape))
    
def is_hermitian(matrix):
    """Determine if a matrix is hermitian."""
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
    return array.__array_wrap__(np.fft.ifftshift(array_hermitian.view(np.ndarray)))

def expand_aperture(ap):
    """Expand an aperture, illuminating neighboring points."""
    nearby = (ap != 0)
    nearby[1:-1,1:-1] |= ap[:-2,:-2]
    nearby[1:-1,1:-1] |= ap[2:,2:]
    nearby[1:-1,1:-1] |= ap[:-2,2:]
    nearby[1:-1,1:-1] |= ap[2:,:-2]
    return nearby
    
def create_complex_HDU(data, kind=None, header=None):
    """Create a FITS HDU which contains complex-valued data.
    
    Complex valued data is stored by adding a new dimension, with shape 2 which
    contains the real and imaginary part of each complex value. The HDU header
    is updated with a few keywords which describe the nature of the data. The
    keyword ``COMPLEX`` is a boolean identifying a complex-valued HDU. The
    keyword ``CAXIS`` is the axis along which the data is split into complex
    and real-valued parts. ``REC`` is an informational keyword which describes
    how to re-assemble the data in python.

    Parameters
    ----------
    data : array-like
        The data to convert to a complex HDU
    kind : HDU class, optional
        The kind of HDU to return, if not provided, defaults to ImageHDU.
    header : FITS Header, optional
        A header to insert into the HDU
    
    Returns
    -------
    HDU : FITS HDU
        The FITS Header-Data Unit with complex data.
    
    """
    from astropy.io import fits
    from astropy.io.fits.hdu.base import _BaseHDU
    cls = fits.ImageHDU if kind is None else kind
    if not issubclass(cls, _BaseHDU):
        raise ValueError("HDUs must be valid FITS HDU kinds, got {0!r}.".format(cls))
    HDU = cls(np.array([data.real, data.imag]), header=header)
    HDU.header["COMPLEX"] = True, "Data is complex"
    HDU.header["CAXIS"] = data.ndim + 1, "Complex data axis"
    HDU.header["REC"] = ("data[0]+1j*data[1]", "To reconstruct complex data in python")
    return HDU
    
def read_complex_HDU(HDU, force=False):
    """Read data from a complex HDU.
    
    Loads data from an HDU, and converts split complex data (such as data
    written by :func:`create_complex_HDU`) into a numpy complex array. Uses
    FITS keywords to identify complex HDUs. ``COMPLEX`` should be a boolean
    which identifies a complex HDU. ``CAXIS`` should identify the FITS axis
    along which the data is split between real and complex values. Headers
    which don't appear to have complex data (or have ``COMPLEX`` set to
    ``False``) will just return the data from the HDU.
    
    When ``CAXIS`` is missing, the complex axis is assumed to be the last FITS
    axis (alternatively the first numpy axis.)
    
    Parameters
    ----------
    HDU : a FITS HDU
        The HDU which might have complex data.
    force : bool, optional
        Whether to force the HDU to appear to have complex data.
    
    """
    if HDU.header.get('COMPLEX', False) or force:
        ndim = int(HDU.header["NAXIS"])
        caxis = int(HDU.header.get("CAXIS", ndim - 1))
        caxis = ndim - caxis
        cslice = [ slice(None, None) for i in range(ndim)]
        cslice[caxis] = 0
        creal = tuple(cslice)
        cslice[caxis] = 1
        cimag = tuple(cslice)
        return HDU.data[creal] + 1j * HDU.data[cimag]
    else:
        return HDU.data