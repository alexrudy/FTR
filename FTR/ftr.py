# -*- coding: utf-8 -*-
"""
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
import weakref

# Scientific Python Imports
import numpy as np
try:
    import scipy.fftpack as fftpack
except ImportError:
    import numpy.fft as fftpack

from astropy.utils import lazyproperty

# Local imports
from .base import Reconstructor
from .io import IOBase
from .libftr._ftr import CFTRBase
from .utils import (complexmp, ignoredivide, remove_piston, remove_tiptilt, 
    fftgrid, shapegrid, shapestr, apply_tiptilt, create_complex_HDU,
    read_complex_HDU)

__all__ = ['FTRFilter', 'FourierTransformReconstructor', 'FastFTReconstructor', 'mod_hud_filter', 'fried_filter', 'ideal_filter', 'inplace_filter', 'hud_filter']

class FTRFilter(collections.namedtuple("_FTRFilter", ["gx", "gy", "name"]), IOBase):
    """
    FTRFilter(gx,gy,name)

    Tuple collection of filter values.
    """
    
    def __to_hdf5__(self, filename, **kwargs):
        """Write an HDF5 file."""
        import h5py
        with h5py.File(filename, kwargs.pop('mode', 'w')) as file:
            group = file.create_group("FTR Filter")
            group.attrs['name'] = self.name
            for g in ('gx', 'gy'):
                dataset = group.create_dataset(g, data=getattr(self, g), **kwargs)
                dataset.attrs['description'] = 'Filter {0:s} component'.format(g.lstrip('g'))
        
    @classmethod
    def __from_hdf5__(cls, filename, **kwargs):
        """Load a filter from an HDF5 file."""
        import h5py
        clsargs = {}
        with h5py.File(filename, kwargs.pop('mode', 'r')) as file:
            group = file['FTR Filter']
            clsargs['name'] = group.attrs['name']
            for name, dataset in group.items():
                clsargs[name] = dataset[...]
        return cls(**clsargs)
    
    def __to_fits__(self, filename, **kwargs):
        """Write the FITS file."""
        from astropy.io import fits
        HDUs = []
        for g in ('gx', 'gy'):
            gHDU = create_complex_HDU(getattr(self,g))
            gHDU.name = g.upper()
            gHDU.header['GAXIS'] = g.lstrip('g').upper()
            gHDU.header['DATA'] = (g.upper(), "Filter {0:s} component".format(g.lstrip('g')))
            gHDU.header['FNAME'] = (self.name, "Filter Name")
            HDUs.append(gHDU)
        fits.HDUList([fits.PrimaryHDU(HDUs[0].data, HDUs[0].header)]+HDUs[1:]).writeto(filename, **kwargs)
    
    @classmethod
    def __from_fits__(cls, filename, **kwargs):
        """Read from a FITS file."""
        from astropy.io import fits
        clsargs = {}
        with fits.open(filename, **kwargs) as HDUs:
            for HDU in HDUs:
                gdata = read_complex_HDU(HDU)
                gaxis = HDU.header.get('GAXIS')
                clsargs['g{0:s}'.format(gaxis.lower())] = gdata
                clsargs['name'] = HDU.header["FNAME"]
        return cls(**clsargs)

class FourierTransformReconstructor(Reconstructor):
    """A reconstructor which uses the fourier transform to turn slopes into an
    estiate of the wavefront phase.
    
    Parameters
    ----------
    ap: array_like
        The aperture of valid measurement points, which also defines the
        reconstructor shape used to generate filters.
    filter: string_like, optional
        The filter name to use. If not provided, it is expected that the user
        will initialize :attr:`gx` and :attr:`gy` themselves.
    manage_tt: bool, optional
        Remove tip and tilt from slopes before reconstruction, and re-apply
        them after reconstruction. (default is to use ``suppress_tt``)
    suppress_tt: bool, optional
        Remove tip and tilt from slopes, and don't re-apply after
        reconstruction. (default is False)
        
    
    Notes
    -----
    
    The Fourier Transform Reconstructor implements the reconstruction scheme
    described in chapter 2 of Lisa Poyneer's dissertation. [1]_ The
    implementation of the Fourier Transform Reconstructor depends on a pair of
    spatial filters (x and y) which relate the specific geometry of the
    wavefront sensor to the phase estimation points. Common spatial filters are
    the "modified Hudgins" filter, ``mod_hud``, and the Fried geometry
    ``fried``. Additional named filters can be registered with this class using
    :meth:`~FourierTransformReconstructor.register`, or an entirely custom
    filter can be used by modifying the :attr:`gx` and :attr:`gy` attributes.
    
    
    References
    ----------
    
    .. [1] Poyneer, L. A. Signal processing for high-precision wavefront
       control in adaptive optics. (Thesis (Ph.D.) - University of California,
       2007).
    
    Examples
    --------
    
    Creating a generic reconstructor will result in an uninitialized filter::
        
        >>> import numpy as np
        >>> aperture = np.ones((10,10))
        >>> recon = FourierTransformReconstructor(aperture)
        >>> recon
        <FourierTransformReconstructor (10x10) filter='Unknown'>
        
    
    You can create a reconstructor with a named filter::
        
        >>> import numpy as np
        >>> aperture = np.ones((10,10))
        >>> recon = FourierTransformReconstructor(aperture, "fried")
        >>> recon
        <FourierTransformReconstructor (10x10) filter='fried'>
        >>> ys, xs = np.meshgrid(np.arange(10), np.arange(10))
        >>> recon(xs, ys)
        array([...])
    
    Reconstructor filters can be changed dynamically::
        
        >>> import numpy as np
        >>> aperture = np.ones((10,10))
        >>> recon = FourierTransformReconstructor(aperture, "fried")
        >>> recon
        <FourierTransformReconstructor (10x10) filter='fried'>
        >>> recon.use("mod_hud")
        >>> recon
        <FourierTransformReconstructor (10x10) filter='mod_hud'>
        >>> recon.name = "hud"
        >>> recon
        <FourierTransformReconstructor (10x10) filter='hud'>
    
    You can set filter components individually::
        
        >>> import numpy as np
        >>> aperture = np.ones((10,10))
        >>> recon = FourierTransformReconstructor(aperture, "fried")
        >>> recon
        <FourierTransformReconstructor (10x10) filter='fried'>
        >>> recon.gx = np.ones((10,10), dtype=np.complex) + 1j * np.sin(np.arange(10) / 10 * np.pi)
        >>> recon.name = "fried custom-x"
        >>> recon
        <FourierTransformReconstructor (10x10) filter='fried custom-x'>
        
    
    """
    
    _gx = None
    _gy = None
    _n = 0
    _dzero = None
    _filtername = "UNDEFINED"
    _denominator = None
    _suppress_tt = None
    _manage_tt = None
    _filters = None
    
    def __repr__(self):
        """Represent this object."""
        return "<{0} {1} filter='{2}'{3}>".format(self.__class__.__name__, shapestr(self.shape), self.name, "" if self.tt_mode == "unmanaged" else " tt='{0}'".format(self.tt_mode))
    
    def __init__(self, ap, filter=None, manage_tt=None, suppress_tt=False):
        self._post_filters = []
        super(FourierTransformReconstructor, self).__init__()
        ap = np.asarray(ap, dtype=np.bool)
        self._shape = ap.shape
        self.ap = ap
        
        self._filtername = "Unknown"
        if filter is not None:
            self.use(six.text_type(filter))
        
        self.suppress_tt = bool(suppress_tt)
        self.manage_tt = manage_tt
        self.y, self.x = shapegrid(self.shape)
        
    @property
    def suppress_tt(self):
        """Tip-tilt suppression flag."""
        return self._suppress_tt
        
    @suppress_tt.setter
    def suppress_tt(self, suppress_tt):
        """Set the tip-tilt suppression flag."""
        if bool(suppress_tt) and self._manage_tt is False:
            raise ValueError("{0:s}: Can't specify manage_tt=False and "
                "suppress_tt=True, as tip/tilt will not be suppressed when "
                "tip-tilt isn't removed from the slopes initially.".format(
                self.__class__.__name__))
        self._suppress_tt = bool(suppress_tt)
    
    @property
    def manage_tt(self):
        """Tip-tilt management flag."""
        if self._manage_tt is None:
            return self._suppress_tt
        return self._manage_tt
        
    @manage_tt.setter
    def manage_tt(self, manage_tt):
        """Set the tip-tilt management flag."""
        if bool(manage_tt) is False and manage_tt is not None and self._suppress_tt:
            raise ValueError("{0:s}: Can't specify manage_tt=False and "
                "suppress_tt=True, as tip/tilt will not be suppressed when "
                "tip-tilt isn't removed from the slopes initially.".format(
                self.__class__.__name__))
        self._manage_tt = bool(manage_tt) if manage_tt is not None else manage_tt
    
    
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
        """The filter components, as a :class:`FTRFilter` tuple."""
        return FTRFilter(self.gx, self.gy, self.name)
        
    @property
    def shape(self):
        """Shape of the reconstructed grid."""
        return self._shape
        
    @property
    def gx(self):
        """The x spatial filter.
        
        :math:`G_{wx}` in the equation implemented in :meth:`apply_filter`.
        """
        return self._gx
        
    @gx.setter
    def gx(self, gx):
        """Set the x filter"""
        self._gx = self._validate_filter(gx)
        self._denominator = None
        
    @property
    def gy(self):
        """The y filter.
        
        :math:`G_{wy}` in the equation implemented in :meth:`apply_filter`.
        """
        return self._gy
        
    @gy.setter
    def gy(self, gy):
        """Set and validate the y filter"""
        self._gy = self._validate_filter(gy)
        self._denominator = None
        
    @property
    def ap(self):
        """The aperture, a boolean numpy array.
        
        The aperture is used to compute tip and tilt removal from the Fourier Transform Reconstructor.
        """
        return self._ap
        
    @ap.setter
    def ap(self, value):
        """Validate aperture."""
        ap = np.asarray(value) != 0
        if ap.shape != self.shape:
            raise ValueError("Aperture should be same shape as input data. Found {0!r}, expected {1!r}.".format(ap.shape, self.shape))
        if not ap.any():
            raise ValueError("Aperture must be illuminated somewhere!")
        self._ap = ap
        
    @property
    def tt_mode(self):
        """Tip/tilt management mode.
        
        Returns a string representation of the tip/tilt management mode, useful
        for the representation of this object.
        """
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
        """Filter denominator
        
        This term normalizes for the magnitude of the individual spatial
        filters. It is recomputed whenever the filters change, but it is
        otherwise computed only once for each set of filters.
        
        The equation used is
        
        .. math::
            
            |G_{wx}|^{2} + |G_{wy}|^2
        
        
        """
        if self._denominator is not None:
            return self._denominator
        self._denominator = (np.abs(self.gx)**2.0 + np.abs(self.gy)**2.0)
        self._denominator[(self._denominator == 0.0)] = 1.0 #Fix non-hermitian parts.
        return self._denominator
        
    def add_post_filter(self, filter):
        """Add a filter to be applied to the phase."""
        self._post_filters.append(filter)
        
    def apply_filters(self, est_ft):
        """Apply phase filters."""
        for filt in self._post_filters:
            est_ft = filt(est_ft)
        return est_ft
        
    def apply_filter(self, xs_ft, ys_ft):
        r"""Apply the filter to the FFT'd values.
        
        Parameters
        ----------
        xs_ft : array_like
            The fourier transform of the x slopes
        ys_ft : array_like
            The fourier transform of the y slopes
            
        Returns
        -------
        est_ft : array_like
            The fourier transform of the phase estimate.
        
        Notes
        -----
        
        This implements the equation
        
        .. math::
            
            \hat{\Phi} = \frac{G_{wx}^{*} X + 
            G_{wy}^{*} Y}{|G_{wx}|^{2} + |G_{wy}|^2}
        
        """
        return ((np.conj(self.gx) * xs_ft + np.conj(self.gy) * ys_ft)
                / self.denominator)
        
    def reconstruct(self, xs, ys, manage_tt=False, suppress_tt=False):
        """Use the Fourier transform and spatial filters to reconstruct an
        estimate of the phase.
        
        Instead of using this method directly, call the instnace itself
        to ensure that settings are correctly obeyed.
        
        Parameters
        ----------
        xs : array_like
            The x slopes
        ys : array_like
            The y slopes
        manage_tt : bool
            Whether to remove the tip/tilt from the slopes before reconstruction
        suppress_tt : bool
            If set, do not re-apply the tip tilt after reconstruction.
        
        Returns
        -------
        estimate : array_like
            An estimate of the phase across all the points
            where x and y slopes were measured.
        
        Notes
        -----
        This method serves as the implementation for :meth:`__call__`.
        
        """
        if manage_tt:
            xs, xt = remove_piston(self.ap, xs)
            ys, yt = remove_piston(self.ap, ys)
        
        xs_ft = fftpack.fftn(xs)
        ys_ft = fftpack.fftn(ys)
        
        # There is a factor of two here. This is applied because the normalization of the FFT in IDL is different from the one used here. See FFTNormalization.rst
        est_ft = self.apply_filter(xs_ft, ys_ft) * 2.0
        est_ft = self.apply_filters(est_ft)
        
        estimate = np.real(fftpack.ifftn(est_ft))
        
        if manage_tt and not suppress_tt:
            estimate = apply_tiptilt(self.ap, estimate, xt, yt)
        
        return estimate
        
    def __call__(self, xs, ys):
        """Perform the reconstruction.
        
        Parameters
        ----------
        xs : array_like
            The x slopes
        ys : array_like
            The y slopes
        
        Returns
        -------
        estimate : array_like
            An estimate of the phase across all the points
            where x and y slopes were measured.
        """
        return self.reconstruct(xs, ys, 
            manage_tt=self.manage_tt, suppress_tt=self.suppress_tt)
        
    def invert(self, estimate):
        """Invert the estimate to produce slopes.
        
        Parameters
        ----------
        estimate : array_like
            Phase estimate to invert.
        
        Returns
        -------
        xs : array_like
            Estimate of the x slopes.
        ys : array_like
            Estimate of the y slopes.
        
        
        """
        if self.manage_tt:
            estimate, ttx, tty = remove_tiptilt(self.ap, estimate)
        
        est_ft = fftpack.fftn(estimate) / 2.0
        
        xs_ft = self.gx * est_ft
        ys_ft = self.gy * est_ft
        
        xs = np.real(fftpack.ifftn(xs_ft))
        ys = np.real(fftpack.ifftn(ys_ft))
        
        if self.manage_tt and not self.suppress_tt:
            xs += ttx
            ys += tty
        
        return (xs, ys)
        
    _REGISTRY = {}
        
    @classmethod
    def register(cls, name, filter=None):
        """Register a filter generating function.
        
        This classmethod can be used as a decorator.
        
        Parameters
        ----------
        name : str
            The name of the filter to register
        filter : callable
            A callable which will return an :class:`FTRFilter`-compatible
            tuple when provided with the desired shape of the filter.
            
        Notes
        -----
        
        To use this method as a decorator, simply provide the filter name::
            
            @FourierTransformReconstructor.register("myfilter")
            def my_filter_function(shape):
                ...
                return FTRFilter(gx, gy, "myfilter")
            
        
        """
        
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
        """Use a particular spatial filter for reconstruction.
        
        Parameters
        ----------
        filter : str
            The filter name, as it was registered with :meth:`register`
        
        Examples
        --------
        
        You can create a reconstructor with one filter::
            
            >>> import numpy as np
            >>> aperture = np.ones((10,10))
            >>> recon = FourierTransformReconstructor(aperture, "fried")
            >>> recon
            <FourierTransformReconstructor (10x10) filter='fried'>
            >>> recon.use("mod_hud")
            >>> recon
            <FourierTransformReconstructor (10x10) filter='mod_hud'>
            
        
        """
        self.gx, self.gy, self._filtername = self.get(filter, self.shape)
        
    @classmethod
    def get(cls, filter, shape):
        """Get a filter tuple by name from the filter registry.
        
        Parameters
        ----------
        filter : str
            The filter name, as it was registered with :meth:`register`
        shape : tuple of ints
            The shape of the desired filter.
        
        Returns
        -------
        gx : array_like
            The x spatial filter
        gy : array_like
            The y spatial filter
        name : str
            The filter name
        
        """
        return cls._REGISTRY[filter](shape)
        
    @classmethod
    def filters(cls):
        """Return the list of registered filter names."""
        return cls._REGISTRY.keys()

class FastFTReconstructor(CFTRBase, FourierTransformReconstructor):
    """A fourier transform reconstructor which implements the reconstruct
    method using :ref:`libftr`.
    
    Parameters
    ----------
    ap: array_like
        The aperture of valid measurement points, which also defines the
        reconstructor shape used to generate filters.
    filter: string_like, optional
        The filter name to use. If not provided, it is expected that the user
        will initialize :attr:`gx` and :attr:`gy` themselves.
    manage_tt: bool, optional
        Remove tip and tilt from slopes before reconstruction, and re-apply
        them after reconstruction. (default is to use ``suppress_tt``)
    suppress_tt: bool, optional
        Remove tip and tilt from slopes, and don't re-apply after
        reconstruction. (default is False)
        
    
    Notes
    -----
    The Fourier transforms used in this reconstructor are handled by FFTW <http://fftw.org>, the Fastest Fourier Transform in the West. FFTW gains much of its speed by operating in-place, and performing as many optimization tricks as it knows how. When benchmarked against the pure-python implementaion above, this implementation should be almost an order of magnitude faster.
    
    See Also
    --------
    :class:`FourierTransformReconstructor`
    
    """
    
    _denominator = None
    
    def reconstruct(self, xs, ys):
        """Reconstruct slopes to phase.
        
        Parameters
        ----------
        xs : array_like
            The x slopes
        ys : array_like
            The y slopes
            
        Returns
        -------
        estimate : array_like
            An estimate of the phase from the slopes.
        
        """
        if self.manage_tt:
            xs, xt = remove_piston(self.ap, xs)
            ys, yt = remove_piston(self.ap, ys)
        
        estimate = super(FastFTReconstructor, self).reconstruct(xs, ys)
        
        if self.manage_tt and not self.suppress_tt:
            estimate = apply_tiptilt(self.ap, estimate, xt, yt)
        
        return estimate
    


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
    # gx[:,nx//2] = 0.0
    # gy[ny//2,:] = 0.0
    
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

    gx[ny//2,:] = 0.0
    gy[:,nx//2] = 0.0
    
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

# @FourierTransformReconstructor.register("ideal")
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
    