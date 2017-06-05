"""
Cython wrappers for the C implementation of the Fourier Transform Reconstructor.
"""

import numpy as np
cimport numpy as np

from ..utils import shapestr

np.import_array()

cdef class CLQGBase:
    
    def __cinit__(self, gains, alphas, hp_coefficients):
        
        # Perform the same initialization checks that python would have done
        # for all of the arrays. These may occur twice (once in __init__),
        # but since this happens first, and could cause problems in the c
        # library, we duplicated it here.
        from ..lqg import _validate_lqg_arguments
        gains, alphas, hp_coefficients = _validate_lqg_arguments(gains, alphas, hp_coefficients)
        
        self._shape = gains.shape[1:]
        self._nlayers = gains.shape[0]
        
        # Now we can set up the attributes, retaining a reference to the filter
        # components so that we own the underlying memory.
        self._gains = gains
        self._alphas = alphas
        self._hp_coefficients = hp_coefficients
        
        # Finally, generate the filter itself. It is an opaque pointer.
        self._filter = lqg_new_filter(self._nlayers, self._shape[0], self._shape[1], 
            <complex *>np.PyArray_DATA(self._gains), <complex *>np.PyArray_DATA(self._alphas), 
            <complex *>np.PyArray_DATA(self._hp_coefficients))
        
    def __dealloc__(self):
        lqg_destroy(self._filter)
        
    def reset(self):
        lqg_reset(self._filter)
        
    def apply_filter(self, est_ft):
        lqg_apply_filter(self._filter, <complex *>np.PyArray_DATA(est_ft))
        return est_ft
        
    def __call__(self, est_ft):
        est_ft = np.asanyarray(est_ft, dtype=np.complex).copy()
        return self.apply_filter(est_ft)
        
    def _preload(self, est_ft):
        """Preload a C-buffer for data processing."""
        self._est_ft_internal_py = np.asanyarray(est_ft, dtype=np.complex).copy()
        self._est_ft_internal_cy = <complex *>np.PyArray_DATA(self._est_ft_internal_py)
    
    def _execute(self):
        """Execute wiht a pre-loaded C buffer."""
        lqg_apply_filter(self._filter, self._est_ft_internal_cy)
    
    