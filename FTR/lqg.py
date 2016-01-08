# -*- coding: utf-8 -*-
"""
Implementation of the LQG reconstructor in python.
"""

import numpy as np

from .ftr import Filter
from .io import IOBase
from .utils import shapestr

__all__ = ['LQGFilter']

class LQGFilter(Filter, IOBase):
    """An LQG Filter container."""
    def __init__(self, gains, alphas, hp_coeffs):
        super(LQGFilter, self).__init__()
        gains = np.asanyarray(gains, dtype=np.complex)
        if not gains.ndim == 3:
            raise ValueError("Gains must have 3 dimensions. gains.ndim={0:d} gains.shape={1!r}".format(gains.ndim, gains.shape))
        self._nlayers = gains.shape[0]
        self._shape = gains.shape[1:]
        self._gains = gains
        self._alphas = np.asanyarray(alphas, dtype=np.complex)
        self._hp_coeffs = np.asanyarray(hp_coeffs, dtype=np.complex)
        
        try:
            assert self._gains.shape == (self._nlayers,) + self._shape, "Gains"
            assert self._alphas.shape == (self._nlayers,) + self._shape, "Alphas"
        except AssertionError:
            raise ValueError("Shape mismatch for {0:s}".format(str(e)))
        
        self.reset()
        
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
        
    def __repr__(self):
        """A string representation for this filter."""
        return "<{0} nlayers={1:d} shape={2:s}>".format(self.__class__.__name__, self.nlayers, shapestr(self.shape))
        
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
        
    def __to_hdf5__(self, file, **kwargs):
        """Write to a HDF5 file."""
        import h5py
        with h5py.File(file, kwargs.pop('mode', 'w')) as file:
            group = file.create_group("LQG Filter")
            for name in [ "gains", "alphas", "highpass_coefficients" ]:
                group.create_dataset(name, data=getattr(self, name), **kwargs)
            
        
    @classmethod
    def __from_hdf5__(cls, file, **kwargs):
        """Read an LQG filter from an HDF5 file."""
        import h5py
        with h5py.File(file, kwargs.pop('mode', 'r')) as file:
            group = file['LQG Filter']
            args = [group[name] for name in [ "gains", "alphas", "highpass_coefficients" ]]
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
        
        HDUs[0].header['NLAYERS'] = self.nlayers, "Number of layers"
        HDUs[0].header['SAXES'] = len(self.shape),
        for i, axis_length in enumerate(self.shape):
            HDUs[0].header['SAXIS{0:d}'.format(i)] = axis_length
        
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
        """Generate an integral controller"""
        if not isinstance(shape, tuple):
            shape = tuple(shape)
        gains = np.float(gain) * np.ones((1,) + shape, dtype=np.complex)
        alphas = np.float(alpha) * np.ones((1,) + shape, dtype=np.complex)
        hp_coeffs = np.zeros(shape, dtype=np.complex)
        return cls(gains, alphas, hp_coeffs)
    

def create_complex_HDU(data):
    """Create a FITS HDU which contains complex-valued data"""
    from astropy.io import fits
    HDU = fits.ImageHDU(np.array([data.real, data.imag]))
    HDU.header["COMPLEX"] = True, "Data is complex"
    HDU.header["CAXIS"] = data.ndim + 1, "Complex data axis"
    HDU.header["REC"] = ("data[0]+1j*data[1]", "To reconstruct complex data")
    return HDU
    
def read_complex_HDU(HDU, force=False):
    """Read data from a complex HDU."""
    from astropy.io import fits
    if HDU.header['COMPLEX'] or force:
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
