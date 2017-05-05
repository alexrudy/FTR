"""
Cython wrappers for the C implementation of the Fourier Transform Reconstructor.
"""

import numpy as np
from ..utils import shapestr
np.import_array()

cdef class CLQGBase:
    
    def __cinit__(self, gains, alphas, hp_coefficients):
        
        # Perform the same initialization checks that python would have done
        # for all of the arrays. These may occur twice (once in __init__),
        # but since this happens first, and could cause problems in the c
        # library, we duplicated it here.
        gains = np.asanyarray(gains, dtype=np.complex)
        alphas = np.asanyarray(alphas, dtype=np.complex)
        hp_coefficients = np.asanyarray(hp_coefficients, dtype=np.complex)
        
        self._shape = gains.shape[1:]
        self._nlayer = gains.shape[0]
        
        if gains.ndim != 3:
            raise ValueError("Gains must have ndim=3, got ndim={0:d}".format(gains.ndim))
        
        try:
            assert gains.shape == (self._nlayer,) + self._shape, "Gains"
            assert alphas.shape == (self._nlayer,) + self._shape, "Alphas"
            assert hp_coefficients.shape == self._shape, "High-pass coefficients"
        except AssertionError as e:
            raise ValueError("Shape mismatch for {0:s}".format(str(e)))
            
        # Now we can set up the attributes, retaining a reference to the filter
        # components so that we own the underlying memory.
        self._gains = gains
        self._alphas = alphas
        self._hp_coefficients = hp_coefficients
        
        # Finally, generate the filter itself. It is an opaque pointer.
        self._filter = lqg_new_filter(self._nlayer, self._shape[0], self._shape[1], 
            <complex *>np.PyArray_DATA(self._gains), <complex *>np.PyArray_DATA(self._alphas), 
            <complex *>np.PyArray_DATA(self._hp_coefficients))
        
    def __dealloc__(self):
        lqg_destroy(self._filter)
        
    def reset(self):
        lqg_reset(self._filter)
        
    def __call__(self, estimate):
        est_ft = np.asanyarray(estimate, dtype=np.complex)
        lqg_apply_filter(self._filter, <complex *>np.PyArray_DATA(est_ft))
        return est_ft
        
    
        
    