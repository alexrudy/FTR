# -*- coding: utf-8 -*-
"""
:mod:`ftr` – The Fourier Transfomr reconstructor
================================================

The fourier transform reconstructor converts slopes (x and y slope grids) to
phase values. The reconstruction works using the fourier transform of the x
and y slopes, and applying a filter, which accounts for the way in which those
slopes map to phase.

The concept for the reconstructor, and the filters documented here, are taken
from Lisa Poyneer's 2007 dissertaiton.

"""

from __future__ import (absolute_import, unicode_literals, division,
                        print_function)

# Python Imports
import abc
import six
import collections

# Scientific Python Imports
import numpy as np
try:
    import scipy.fftpack as fftpack
except ImportError:
    import numpy.fft as fftpack

from astropy.utils import lazyproperty

# Local imports
from .base import Reconstructor
from .utils import (complexmp, ignoredivide, remove_piston, fftgrid, 
    shapegrid, shapestr)

__all__ = ['FTRFilter', 'FourierTransformReconstructor', 'mod_hud_filter', 'fried_filter', 'ideal_filter']

FTRFilter = collections.namedtuple("FTRFilter", ["gx", "gy", "name"])

class FourierTransformReconstructor(Reconstructor):
    """The Fourier Transform Reconstructior uses the fourier domain relationship
    between the gradient and phase.
    
    Parameters
    ----------
    shape: tuple of int
        The size of the reconstruction grid
    ap: array_like, (n x n)
        The aperture of valid measurement points
    filter: string_like, optional
        The filter name to use. If not provided, it is expected that the user
        will initialize :attr:`gx` and :attr:`gy` themselves.
    manage_tt: bool, optional
        Remove tip and tilt from slopes before reconstruction, and re-apply them after reconstruction. (default is to use ``suppress_tt``)
    suppress_tt: bool, optional
        Remove tip and tilt from slopes, and don't re-apply after reconstruction. (default is False)
    
    
    """
    
    _gx = None
    _gy = None
    _n = 0
    _dzero = None
    _filtername = "UNDEFINED"
    
    def __repr__(self):
        """Represent this object."""
        return "<{0} {1} filter='{2}'{3}>".format(self.__class__.__name__, shapestr(self.shape), self.name, "" if self.tt_mode is "unmanaged" else " tt='{0}'".format(self.tt_mode))
    
    def __init__(self, shape, ap, filter=None, manage_tt=None, suppress_tt=False):
        super(FourierTransformReconstructor, self).__init__()
        self._shape = tuple([np.int(s) for s in shape])
        self._filtername = "Unknown"
        if filter is not None:
            self.use(six.text_type(filter))
        self.ap = ap
        
        self.suppress_tt = bool(suppress_tt)
        if manage_tt is None:
            self.manage_tt = bool(suppress_tt)
            self.suppress_tt = bool(suppress_tt)
        else:
            self.manage_tt = bool(manage_tt)
            if suppress_tt and not self.manage_tt:
                raise TypeError("{0:s}: Can't specify manage_tt=False and "
                "suppress_tt=True, as tip/tilt will not be suppressed when "
                "tip-tilt isn't removed from the slopes initially.")
            self.suppress_tt = False
        self.y, self.x = shapegrid(self.shape)
        
    @property
    def name(self):
        """Filter name"""
        return self._filtername
        
    @name.setter
    def name(self, value):
        """Set the filter name."""
        if value in self._REGISTRY:
            self.use(value)
        else:
            self._filtername = value
        
    @property
    def filter(self):
        """The filter components."""
        return FTRFilter(self.gx, self.gy, self.name)
        
    @property
    def shape(self):
        """Shape of the reconstructed grid."""
        return self._shape
        
    @property
    def gx(self):
        """The x filter"""
        return self._gx
        
    @gx.setter
    def gx(self, gx):
        """Set the x filter"""
        self._gx = self._validate_filter(gx)
        self._denominator = None
        
    @property
    def gy(self):
        """The y filter"""
        return self._gy
        
    @gy.setter
    def gy(self, gy):
        """Set and validate the y filter"""
        self._gy = self._validate_filter(gy)
        self._denominator = None
        
    @property
    def ap(self):
        """The aperture"""
        return self._ap
        
    @ap.setter
    def ap(self, value):
        """Validate aperture."""
        ap = np.asarray(value).astype(np.bool)
        if ap.shape != self.shape:
            raise ValueError("Aperture should be same shape as input data. Found {0!r}, expected {1!r}.".format(ap.shape, self.shape))
        if not ap.any():
            raise ValueError("Aperture must be illuminated somewhere!")
        self._ap = ap
        
    @property
    def tt_mode(self):
        """Tip/tilt management mode."""
        if not self.manage_tt:
            return "unmanaged"
        elif self.manage_tt and not self.suppress_tt:
            return "managed"
        elif self.manage_tt and self.suppress_tt:
            return "suppressed"
        
    def _validate_filter(self, _filter):
        """Ensure that a filter is the correct shape.
        
        This method checks the shape and that the filter is all finite.
        
        :param _filter: The filter array to check.
        :returns: The filter, correctly typed and checked for consistency.
        """
        gf = np.asarray(_filter).astype(np.complex)
        if gf.shape != self.shape:
            raise ValueError("Filter should be same shape as input data. Found {0!r}, expected {1!r}.".format(gf.shape, self.shape))
        if not np.isfinite(gf).all():
            raise ValueError("Filter must be finite at all points!")
        return gf
        
    @property
    def denominator(self):
        """Filter denominator"""
        if self._denominator is not None:
            return self._denominator
        self._denominator = np.abs(self.gx)**2.0 + np.abs(self.gy)**2.0
        self._denominator[(self._denominator == 0.0)] = 1.0 #Fix non-hermetian parts.
        return self._denominator
        
    def apply_filter(self, xs_ft, ys_ft):
        """Apply the filter to the FFT'd values.
        
        :param xs_ft: The x fourier transform
        :param ys_ft: THe y fourier transform
        :returns: The filtered estimate, fourier transformed.
        
        """
        return (np.conj(self.gx) * xs_ft + np.conj(self.gy) * ys_ft)/self.denominator
        
    def reconstruct(self, xs, ys):
        """The reconstruction method"""
        if self.manage_tt:
            xs, xt = remove_piston(self.ap, xs)
            ys, yt = remove_piston(self.ap, ys)
        
        xs_ft = fftpack.fftn(xs)
        ys_ft = fftpack.fftn(ys)
        
        est_ft = self.apply_filter(xs_ft, ys_ft)
        
        estimate = np.real(fftpack.ifftn(est_ft))
        
        if self.manage_tt and not self.suppress_tt:
            return estimate + (self.x * xt) + (self.y * yt)
        
        return estimate
        
    _REGISTRY = {}
        
    @classmethod
    def register(cls, name, filter=None):
        """Register a filter generating function."""
        
        def _register(filterfunc):
            """Filter Function"""
            cls._REGISTRY[name] = filterfunc
            return filterfunc
        
        if six.callable(name):
            filterfunc = name
            cls._REGISTRY[filterfunc.__name__] = filterfunc
            return filterfunc
        elif isinstance(name, six.text_type) and filter is None:
            return _register
        elif isinstance(name, six.text_type) and six.callable(filter):
            return _register(filterfunc)
        else:
            raise TypeError("Filter must be a callable, or a name, and used as a decorator.")
    
    def use(self, filter):
        """Use a particular filter."""
        self.gx, self.gy, self._filtername = self._REGISTRY[filter](self.shape)
        
    @classmethod
    def get(cls, filter, shape):
        """Get a filter by name."""
        return cls._REGISTRY[filter](shape)
        
    @classmethod
    def filters(cls):
        """Return the list of filters available."""
        return cls._REGISTRY.keys()

@FourierTransformReconstructor.register("hud")
def hud_filter(shape):
    """Hudgins shearing interferometer. In this geometry the WFS produces
    measurements which are the differences in phase values between two points.
    This corresponds to wavefront sensors centered between each pair of
    points."""
    fy, fx = fftgrid(shape)
    ny, nx = shape
    
    gx = np.exp(1j * fx) - 1.0
    gy = np.exp(1j * fy) - 1.0
    
    #TODO: Check that this nyquist nulling is correct.
    # gx[:,nx/2] = 0.0
    # gy[ny/2,:] = 0.0
    
    #TODO: Check for other null points in the filter.
    
    return FTRFilter(gx, gy, "hud")

@FourierTransformReconstructor.register("mod_hud")
def mod_hud_filter(shape):
    """The modified hudgins filter is a geomoetry similar to
    a Fried geometry, but where the slope measurements, rather than
    being the difference between two adjacent points, the slope is
    taken to be the real slope at the point in the center of four
    phase measurement points."""
    fy, fx = fftgrid(shape)
    ny, nx = shape
    
    gx = np.exp(1j*fy/2)*(np.exp(1j*fx) - 1)
    gy = np.exp(1j*fx/2)*(np.exp(1j*fy) - 1)

    gx[ny/2,:] = 0.0
    gy[:,nx/2] = 0.0
    
    return FTRFilter(gx, gy, "mod_hud")

@FourierTransformReconstructor.register("fried")
def fried_filter(shape):
    """The fried filter is for a system geometry where the
    slope measruement points are taken to be the difference
    between two adjacent points. As such, the slopes are
    reconstructed to be the phase points 1/2 a subaperture
    away from the measured point.

    In this scheme, the slope measruement in the center of 
    four phase measurement points is taken (in the x-direction)
    to be the average of the x slope between the top measruements
    and the x slope between the bottom measurement."""
    fy, fx = fftgrid(shape)
    ny, nx = shape
    
    gx = (np.exp(1j*fy) + 1)*(np.exp(1j*fx) - 1)
    gy = (np.exp(1j*fx) + 1)*(np.exp(1j*fy) - 1)
    
    gx[ny//2,:] = 0.0
    gy[:,nx//2] = 0.0
    
    return FTRFilter(gx, gy, "fried")

@FourierTransformReconstructor.register("ideal")
def ideal_filter(shape):
    """An Ideal filter represents a phase where the slope
    measurements are taken to be a continuous sampling of
    the phase between phase measurement points. """
    fy, fx = fftgrid(shape)
    ny, nx = shape
    
    with ignoredivide():
        gx = (
                  1/fy  * ( (np.cos(fy) - 1) * np.sin(fx) 
                          + (np.cos(fx) - 1) * np.sin(fy) ) +
                  1j/fy * ( (np.cos(fx) - 1) * (np.cos(fy) - 1) 
                          - np.sin(fx) * np.sin(fy) ) )
        gy = (
                  1/fx  * ( (np.cos(fx) - 1) * np.sin(fy) 
                          + (np.cos(fy) - 1) * np.sin(fx) ) +
                  1j/fx * ( (np.cos(fy) - 1) * (np.cos(fx) - 1) 
                          - np.sin(fy) * np.sin(fx) ) )
      
    # Exclude division by zero!
    gx[0,:] = 0
    gy[:,0] = 0
    
    # the filter is anti-Hermitian here. The real_part takes care
    # of it, but simpler to just zero it out.
    gx[ny//2,:] = 0.0
    gy[:,nx//2] = 0.0
    
    return FTRFilter(gx, gy, "ideal")
    
@FourierTransformReconstructor.register("inplace")
def inplace_filter(shape):
    """This filter modifies the 'fried' filter so that it reconstructs to the in-place positions of the slopes."""
    fy, fx = fftgrid(shape)
    xshift, yshift = (-0.5, -0.5)
    
    gx, gy, name = fried_filter(shape)
    
    gx *= complexmp(1.0, fx * xshift) * complexmp(1.0, fy * yshift)
    gy *= complexmp(1.0, fx * xshift) * complexmp(1.0, fy * yshift)
    
    return FTRFilter(gx, gy, "inplace")
    