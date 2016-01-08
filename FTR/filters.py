# -*- coding: utf-8 -*-
"""
Fourier domain filters for use with FTR.
"""

from __future__ import absolute_import
import numpy as np

from .ftr import Filter
from .io import IOBase
from .utils import create_complex_HDU, read_complex_HDU, shapestr

class StaticFilter(Filter, IOBase):
    """A static-valued filter."""
    def __init__(self, array):
        super(StaticFilter, self).__init__()
        self.array = np.asarray(array, dtype=np.complex)
        
    def apply_filter(self, est_ft):
        """Apply the filter."""
        return est_ft * self.array
        
    @property
    def shape(self):
        """Return the array shape."""
        return self.array.shape
    
    def __repr__(self):
        """Represenetation."""
        return "<{0} {1} at {2}>".format(self.__class__.__name__, shapestr(self.shape), hex(id(self)))
        
    def __to_fits__(self, file, **kwargs):
        """Write to a FITS file."""
        from astropy.io import fits
        kwargs.setdefault('clobber', True)
        
        HDU = create_complex_HDU(self.array)
        HDU.name = "FILTER"
        
        HDU.writeto(file, **kwargs)
        
    @classmethod
    def __from_fits__(cls, file, **kwargs):
        """Read from a FITS file."""
        from astropy.io import fits
        force = kwargs.pop('force', False)
        
        with fits.open(file, **kwargs) as HDUs:
            filter = read_complex_HDU(HDUs['FILTER'])
        return cls(filter)