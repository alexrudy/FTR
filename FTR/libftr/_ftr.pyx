"""
Cython wrappers for the C implementation of the Fourier Transform Reconstructor.
"""

import numpy as np
cimport numpy as np

include "_ftr.pxd"

from ..utils import shapestr

#Need to ensure that we've initialized the numpy c-API.
np.import_array()

cdef class CFTRBase:
    
    cdef ftr_plan _plan
    cdef np.ndarray sx, sy, est
    cdef np.ndarray _gx, _gy
    cdef tuple _shape
    
    def __cinit__(self, ap, *args, **kwargs):
        """
        Initialize the c-allocated variables, including interanl numpy arrays.
        """
        shape = np.asarray(ap).shape
        self._shape = shape
        
        # Set up the arrays here so that we never have to allocate memory again.
        # but use Numpy so that we get reference counting, etc., for free.
        self.sx = np.empty(shape, dtype=np.float)
        self.sy = np.empty(shape, dtype=np.float)
        self.est = np.empty(shape, dtype=np.float)
        
        # We can safely set up the plan now, as we've initialized all
        # of the required vectors.
        self._plan = ftr_plan_reconstructor(shape[0], shape[1], 
            <double *>np.PyArray_DATA(self.sx), <double *>np.PyArray_DATA(self.sy), 
            <double *>np.PyArray_DATA(self.est))
        
        # We track these here independently of the opaque
        # internal struct so that we can return them as properties.
        self._gx = np.empty(shape, dtype=np.complex)
        self._gy = np.empty(shape, dtype=np.complex)
        
        # We still need to call self.update_filters when the values
        # in gx or gy change, so that we can recompute the denominator.
        self.update_filters()
        
    
    def __dealloc__(self):
        """Deallocation should remove the plan, as that is the non-gc'd object."""
        ftr_destroy(self._plan)
    
    def reconstruct(self, np.ndarray[np.float64_t, ndim=2] sx, np.ndarray[np.float64_t, ndim=2] sy):
        """Perform the actual reconstruction, with a single memory copy."""
        self.sx[...] = sx
        self.sy[...] = sy
        ftr_reconstruct(self._plan)
        return self.est
        
    cdef void update_filters(self):
        """Update the filters in the plan given the values for the x and y filter stored in this object.
        
        This function must be called every time gx or gy is altered.
        """
        ftr_set_filter(self._plan, <complex *>np.PyArray_DATA(self._gx), <complex *>np.PyArray_DATA(self._gy))
        
    def __repr__(self):
        return "<CFTRBase {0:s}>".format(shapestr(self.shape))
        
    def __call__(self, sx, sy):
        return self.reconstruct(sx, sy)
    
    property shape:
        def __get__(self):
            return self._shape
    
    property gx:
        def __get__(self):
            return self._gx
        
        def __set__(self, value):
            self._gx[:] = np.asarray(value, dtype=np.complex)
            self.update_filters()
            
    property gy:
        def __get__(self):
            return self._gy

        def __set__(self, value):
            self._gy[:] = np.asarray(value, dtype=np.complex)
            self.update_filters()
            
        
    

cdef class HalfComplexMapping:
    
    cdef ftr_halfcomplex _hcmap
    cdef tuple _shape
    
    def __cinit__(self, shape):
        self._shape = tuple(shape)
        self._hcmap = ftr_halfcomplex_map(self._shape[0], self._shape[1])
        
    def __dealloc__(self):
        ftr_halfcomplex_destroy(self._hcmap)
        
    property full_to_halfcomplex:
        def __get__(self):
            cdef np.npy_intp shape[2]
            shape[0] = <np.npy_intp> self._shape[0]
            shape[1] = <np.npy_intp> self._shape[1]
            return np.PyArray_SimpleNewFromData(1, shape, np.NPY_INT, <void*>self._hcmap.f2hc)
    
    property halfcomplex_to_full:
        def __get__(self):
            cdef np.npy_intp shape[2]
            shape[0] = <np.npy_intp> self._shape[0]
            shape[1] = <np.npy_intp> self._hcmap.nf
            return np.PyArray_SimpleNewFromData(1, shape, np.NPY_INT, <void*>self._hcmap.hc2f)
            
    