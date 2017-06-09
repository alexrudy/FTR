# -*- coding: utf-8 -*-
"""
Implementation of the LQG reconstructor in python.
"""

import numpy as np
import abc

from .base import Filter
from .io import IOBase
from .utils import shapestr, create_complex_HDU, read_complex_HDU
from .libftr.lqg import CLQGBase

__all__ = ['LQGFilter', 'FastLQGFilter']

def _validate_lqg_arguments(gains, alphas, hp_coeffs, dtype=np.complex):
    """Validate LQG arguments"""
    gains = np.asanyarray(gains, dtype=dtype)
    alphas = np.asanyarray(alphas, dtype=dtype)
    hp_coeffs = np.asanyarray(hp_coeffs, dtype=dtype)
    
    if not gains.ndim >= 2:
        raise ValueError("Gains must have more than 1 dimension. gains.ndim={0:d} gains.shape={1!r}".format(gains.ndim, gains.shape))
    nlayers = gains.shape[0]
    shape = gains.shape[1:]
    
    if gains.shape != (nlayers,) + shape:
        raise ValueError("Shape mismatch for {0:s}, got {1!r} expected {2!r}".format('gains', gains.shape, (nlayers,) + shape))
    if alphas.shape != (nlayers,) + shape:
        raise ValueError("Shape mismatch for {0:s}, got {1!r} expected {2!r}".format('alphas', alphas.shape, (nlayers,) + shape))
    if hp_coeffs.shape != shape:
        raise ValueError("Shape mismatch for {0:s}, got {1!r} expected {2!r}".format('hp_coeffs', hp_coeffs.shape, shape))
    return gains, alphas, hp_coeffs
    
class LQGBase(Filter, IOBase):
    """An LQG Filter container which specifies the interface."""
    def __init__(self, gains, alphas, hp_coeffs):
        super(LQGBase, self).__init__()
        self.reset()
        
    @abc.abstractmethod
    def reset(self):
        """Reset the filter."""
        pass
        
    def __repr__(self):
        """A string representation for this filter."""
        return "<{0} nlayers={1:d} shape={2:s}>".format(self.__class__.__name__, self.nlayers, shapestr(self.shape))
    
    @abc.abstractproperty
    def nlayers(self):
        """Number of layers."""
        return self._nlayers
    
    @abc.abstractproperty
    def shape(self):
        """Shape of the filter, without the layer axis."""
        return self._shape
    
    @abc.abstractproperty
    def gains(self):
        """Filter gains."""
        return self._gains
        
    @abc.abstractproperty
    def alphas(self):
        """Filter Alphas"""
        return self._alphas
        
    @abc.abstractproperty
    def highpass_coefficients(self):
        """High-pass coefficients."""
        return self._hp_coeffs
        
    def __to_hdf5__(self, file, **kwargs):
        """Write to a HDF5 file."""
        import h5py
        with h5py.File(file, kwargs.pop('mode', 'w')) as h5file:
            group = h5file.create_group(kwargs.pop("group", "LQG Filter"))
            for name in [ "gains", "alphas", "highpass_coefficients" ]:
                group.create_dataset(name, data=getattr(self, name), **kwargs)
            
        
    @classmethod
    def __from_hdf5__(cls, file, **kwargs):
        """Read an LQG filter from an HDF5 file."""
        import h5py
        with h5py.File(file, kwargs.pop('mode', 'r')) as h5file:
            # Get either the group with the right name, or the first group from the HDF5 file.
            group = h5file.get(kwargs.pop("group", "LQG Filter"), next(iter(h5file.values())))
            args = [group[name][...] for name in [ "gains", "alphas", "highpass_coefficients" ]]
        return cls(*args)
        
    def __to_fits__(self, file, **kwargs):
        """Write to a FITS file."""
        from astropy.io import fits
        
        kwargs.setdefault('clobber', True)
        
        names =    [ "gains", "alphas", "highpass_coefficients" ]
        HDUnames = [ "GAINS", "ALPHAS", "HPCOEFFS"]
        HDUs = [ create_complex_HDU(getattr(self, name)) for name in names ]
        for HDU, name in zip(HDUs, HDUnames):
            HDU.name = name
        
        data_names = ['gains', 'alphas', 'hp coeffs']
        comments =  ['LGQ filter gains', 'LGQ filter alphas', 'High-pass Filter Coefficients']
        for HDU, data, comment in zip(HDUs, data_names, comments):
            HDU.header['DATA'] = (data, comment)
        
        # Apply only to the first two HDUs.
        for HDU in HDUs[0:2]:
            HDU.header['NLAYERS'] = self.nlayers, "Number of layers"
            HDU.header['SAXES'] = len(self.shape), "Shape of the phase grid"
            for i, axis_length in enumerate(self.shape):
                HDU.header['SAXIS{0:d}'.format(i)] = axis_length, "Axis {0} of the phase grid".format(i)
        
        HDUs[0] = fits.PrimaryHDU(HDUs[0].data, HDUs[0].header)
        
        HDUs = fits.HDUList(HDUs)
        HDUs.writeto(file, **kwargs)
        
    @classmethod
    def __from_fits__(cls, file, **kwargs):
        """Read from a FITS file."""
        from astropy.io import fits
        
        force = kwargs.pop('force', False)
        with fits.open(file, **kwargs) as HDUs:
            gains = read_complex_HDU(HDUs[0], force=force)
            alphas = read_complex_HDU(HDUs[1], force=force)
            hp_coeffs = read_complex_HDU(HDUs[2], force=force)
        return cls(gains, alphas, hp_coeffs)
        
    @classmethod
    def generate_integrator(cls, gain, alpha, shape):
        """Generate an integral controller, done by zeroing out the highpass term."""
        if not isinstance(shape, tuple):
            shape = tuple(shape)
        gains = np.float(gain) * np.ones((1,) + shape, dtype=np.complex)
        alphas = np.float(alpha) * np.ones((1,) + shape, dtype=np.complex)
        hp_coeffs = np.zeros(shape, dtype=np.complex)
        return cls(gains, alphas, hp_coeffs)
        
    @classmethod
    def generate_kalman_filter(cls, k1, alpha, shape):
        """Generate a generic Kalman filter."""
        if not isinstance(shape, tuple):
            shape = tuple(shape)
        gains = np.float(alpha * k1) * np.ones((1,) + shape, dtype=np.complex)
        alphas = np.float(alpha) * np.ones((1,) + shape, dtype=np.complex)
        hp_coeffs = np.float(k1) * np.ones(shape, dtype=np.complex)
        return cls(gains, alphas, hp_coeffs)

class LQGFilter(LQGBase):
    """An LQG Filter container."""
    
    def __init__(self, gains, alphas, hp_coeffs):
        gains, alphas, hp_coeffs = _validate_lqg_arguments(gains, alphas, hp_coeffs, dtype=np.complex)
        self._nlayers = gains.shape[0]
        self._shape = gains.shape[1:]
        self._gains = gains
        self._alphas = alphas
        self._hp_coeffs = hp_coeffs
        super(LQGFilter, self).__init__(gains, alphas, hp_coeffs)
        
    
    def reset(self):
        """Reset the filter memory."""
        self._memory_end = np.zeros(self.shape, dtype=np.complex)
        self._memory_past = np.zeros((self.nlayers,) + self.shape, dtype=np.complex)
        
    def apply_filter(self, est_ft):
        """Apply the Kalman Filter."""
        layer_update = est_ft[None,...] * self.gains + self._memory_past * self.alphas
        # Note that we put the [None,...] into est_ft to add an axis along the front. 
        # This new axis corresponds to the number of layers in a given Kalman Filter.
        
        est_ft = layer_update.sum(axis=0) - self._memory_end * self.highpass_coefficients
        # We collapse layer_update along the number of layers axis.
        
        # At the end, we re-insert data into its original holders.
        self._memory_end[...] = est_ft
        self._memory_past[...] = layer_update
        return est_ft
    
    @property
    def nlayers(self):
        """Number of layers."""
        return self._nlayers
    
    @property
    def shape(self):
        """Shape of the filter, without the layer axis."""
        return self._shape
    
    @property
    def gains(self):
        """Filter gains."""
        return self._gains
        
    @property
    def alphas(self):
        """Filter Alphas"""
        return self._alphas
        
    @property
    def highpass_coefficients(self):
        """High-pass coefficients."""
        return self._hp_coeffs

class FastLQGFilter(CLQGBase, LQGBase):
    """A c-based implementation of the LQG filter."""
    
    @property
    def nlayers(self):
        """Number of layers."""
        return self._nlayers
    
    @property
    def shape(self):
        """Shape of the filter, without the layer axis."""
        return self._shape
    
    @property
    def gains(self):
        """Filter gains."""
        return self._gains
        
    @property
    def alphas(self):
        """Filter Alphas"""
        return self._alphas
        
    @property
    def highpass_coefficients(self):
        """High-pass coefficients."""
        return self._hp_coefficients

