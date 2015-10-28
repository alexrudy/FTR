# -*- coding: utf-8 -*-
"""
Tools for managing half-complex data from FFTW.
"""

import numpy as np
cimport numpy as np

include "libftr/_ftr.pxd"

np.import_array()

cdef class HalfComplexMapping:
    """A class for managing the mapping between the half-complex and full Fourier formats used by FFTW.
    
    Parameters
    ----------
    shape : tuple of ints
        The shape of the full-data array to be used. The shape of the half-complex array
    
    """
    
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
            return np.PyArray_SimpleNewFromData(2, shape, np.NPY_INT, <void*>self._hcmap.f2hc)
    
    property halfcomplex_to_full:
        def __get__(self):
            cdef np.npy_intp shape[2]
            shape[0] = <np.npy_intp> self._shape[0]
            shape[1] = <np.npy_intp> self._hcmap.nf
            return np.PyArray_SimpleNewFromData(2, shape, np.NPY_INT, <void*>self._hcmap.hc2f)
            
    property conjugate:
        def __get__(self):
            should_conj = np.ones(self._shape, dtype=np.bool)
            should_conj.flat[self.full_to_halfcomplex] = False
            return should_conj
            
    
    property halfcomplex_shape:
        def __get__(self):
            return (self._shape[0], self._hcmap.nf)
        
    property shape:
        def __get__(self):
            return self._shape
    
    def unpack(self, halfcomplex_data, full_shape):
        halfcomplex_data = np.asarray(halfcomplex_data)
        full_data = np.empty(full_shape, dtype=halfcomplex_data.dtype)
        full_data.flat[:] = np.conj(halfcomplex_data.flat[self.full_to_halfcomplex.flat[:]])
        full_data.flat[self.halfcomplex_to_full] = halfcomplex_data.flat[:]
        return full_data
        
    def pack(self, full_data):
        full_data = np.asarray(full_data)
        if full_data.shape != self.shape:
            raise ValueError("Shape mismatch between mapping {!r} and data {!r}".format(self.shape, full_data.shape))
        halfcomplex_shape = list(full_data.shape)
        halfcomplex_shape[1] = (full_data.shape[1] // 2) + 1
        halfcomplex_data = np.empty(halfcomplex_shape, dtype=full_data.dtype)
        halfcomplex_data.flat[...] = full_data.flat[self.halfcomplex_to_full.flat[:]]
        return halfcomplex_data

def _prepare_arrays(input_data, output_shape):
    """Prepare input and output arrays."""
    input_data  = np.asarray(input_data)
    input_shape = input_data.shape
    output_data = np.empty(output_shape, dtype=input_data.dtype)
    
    # Total number of elements in each type of array.
    nft = np.prod(input_shape)
    nn  = np.prod(output_shape)
    
    # Extract single dimension parameters, where relevant.
    ny = output_shape[0]
    nf = input_shape[-1]
    nx = output_shape[-1]
    if ny != input_shape[0]:
        raise ValueError("Shape mismatch: ({}x ) -> ({}x )".format(input_shape[0], ny))
    if (nx // 2) + 1 != nf:
        raise ValueError("Shape mismatch: ( x{}) -> ( x{}), expected ( x{})".format(nf, nx, (nx // 2) + 1))
    return output_data
    
def unpack_halfcomplex(input_data, output_shape):
    """Unpack a halfcomplex array into a given output shape.
    
    Parameters
    ----------
    input_data : array_like, N x M
        Data to unpack from the half complex format.
    output_shape : tuple, (N, K)
        Output data shape. ``M = (K // 2) + 1``
    
    Returns
    -------
    output_data : array_like, N x K
        Unpacked full data from half complex input.
        
    """
    input_data  = np.asarray(input_data)
    output_data = _prepare_arrays(input_data, output_shape)
    input_shape = input_data.shape
    
    # Extract single dimension parameters, where relevant.
    ny = output_shape[0]
    nf = input_shape[-1]
    nx = output_shape[-1]
    o = 1 if nx % 2 == 1 else 2
    output_data[ :  , 0:nf] = input_data
    output_data[0   ,nf:nx] = np.conj(input_data[0,nf-o:0:-1])
    output_data[1:ny,nf:nx] = np.conj(input_data[ny-1:0:-1,nf-o:0:-1])
    return output_data

def pack_halfcomplex(input_data):
    """Pack halfcomplex data from full data.
    
    Parameters
    ----------
    input_data : array_like, N x K
        Data to unpack from the half complex format.
    
    Returns
    -------
    output_data : array_like, N x M
        Unpacked full data from half complex input. ``M = (K // 2) + 1``
    
    """
    ny = input_data.shape[0]
    nx = input_data.shape[1]
    nf = (nx // 2) + 1
    return input_data[:,:nf].copy()