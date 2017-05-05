"""
Cython wrappers for the C implementation of the Fourier Transform Reconstructor.
"""
cimport numpy as np

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

cdef class CLQGBase:
    
    cdef readonly int _nlayers
    cdef readonly tuple _shape
    cdef lqg_filter _filter
    cdef readonly np.ndarray _aperture, _gains, _alphas, _hp_coefficients

    cdef complex * _est_ft_internal_cy
    cdef np.ndarray _est_ft_internal_py
    
    