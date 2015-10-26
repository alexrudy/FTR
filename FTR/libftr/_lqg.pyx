"""
Cython wrappers for the C implementation of the Fourier Transform Reconstructor.
"""

import numpy as np
cimport numpy as np

from ..utils import shapestr

cdef extern from "complex.h":
    pass
    
cdef extern from "fftw3.h":
    pass

cdef extern from "lqg.h":
    
    cdef struct lqg_filter_s:
        pass
    ctypedef lqg_filter_s *lqg_filter
    
    lqg_filter lqg_new_filter(const int nl, const int ny, const int nx, complex * gains, complex * alphas, complex * hp_coefficients)
    
    void lqg_apply_filter(lqg_filter filter, complex * est_ft)
    void lqg_filter_callback(void * filter, const int ny, const int nx, complex * est_ft)
    void lqg_reset(lqg_filter filter)
    void lqg_destroy(lqg_filter filter)
    
np.import_array()

cdef class CLQGBase:
    
    cdef int _nlayer
    cdef tuple _shape
    cdef lqg_filter _filter
    cdef np.ndarray aperture, gains, alphas, hp_coefficients
    
    def __cinit__(self, nlayer, aperture, gains, alphas, hp_coefficients):
        self._nlayer = nlayer
        aperture = np.asanyarray(aperture, dtype=np.int)
        self._shape = aperture.shape
        self.aperture = aperture
        
        gains = np.asanyarray(gains, dtype=np.complex)
        alphas = np.asanyarray(alphas, dtype=np.complex)
        hp_coefficients = np.asanyarray(hp_coefficients, dtype=np.complex)
        try:
            assert gains.shape == (nlayer,) + self._shape, "Gains"
            assert alphas.shape == (nlayer,) + self._shape, "Alphas"
            assert hp_coefficients.shape == self._shape, "High-pass coefficients"
        except AssertionError as e:
            raise ValueError("Shape mismatch for {0:s}".format(str(e)))
        self.gains = gains
        self.alphas = alphas
        self.hp_coefficients = hp_coefficients
        self._filter = lqg_new_filter(self._nlayer, self._shape[0], self._shape[1], 
            <complex *>np.PyArray_DATA(self.gains), <complex *>np.PyArray_DATA(self.alphas), 
            <complex *>np.PyArray_DATA(self.hp_coefficients))
        
    def __dealloc__(self):
        lqg_destroy(self._filter)
        
    def reset(self):
        lqg_reset(self._filter)
        
    def __call__(self, estimate):
        est_ft = np.asanyarray(estimate, dtype=np.complex)
        lqg_apply_filter(self._filter, <complex *>np.PyArray_DATA(est_ft))
        return est_ft
        
    
        
    