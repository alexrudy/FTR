"""
Cython wrappers for the C implementation of the Fourier Transform Reconstructor.
"""

import numpy as np
cimport numpy as np

cdef extern from "complex.h":
    pass
    
cdef extern from "fftw3.h":
    pass

# Import declarations from the FTR header files.
cdef extern from "ftr.h":
    
    # Leave the FTR plan as opaque to the user for now
    # We might have to change this in the future to support
    # things like fourier-space gains etc.
    cdef struct FTR_plan:
        pass
    ctypedef FTR_plan *ftr_plan
    
    # These four methods create the functional interface to FTR in c.
    ftr_plan ftr_plan_reconstructor(int nx, int ny, double *sx, double *sy, double *est)
    void ftr_set_filter(ftr_plan recon, complex *gx, complex *gy)
    void ftr_reconstruct(ftr_plan recon)
    void ftr_destroy(ftr_plan recon)
    
    # These are the half-complex mapping tools.
    cdef struct ftr_halfcomplex_s:
        int ny, nx, nf
        int *f2hc
        int *hc2f
    ctypedef ftr_halfcomplex_s *ftr_halfcomplex
    
    ftr_halfcomplex ftr_halfcomplex_map(const int ny, const int nx)
    void ftr_halfcomplex_destroy(ftr_halfcomplex ftr_hc)