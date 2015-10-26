"""
Cython wrappers for the slope management tools.
"""

import numpy as np
cimport numpy as np

from ..utils import shapestr

# Import declarations from the SlopeManage header files.
cdef extern from "slopemanage.h":
    
    cdef struct SM_plan:
        pass
    ctypedef SM_plan sm_plan
    
    sm_plan slope_management_plan(int nx, int ny, int *ap)
    void slope_management_execute(sm_plan plan, double * sx, double * sy)
    void slope_management(int nx, int ny, int *ap, double * sx, double * sy)
    void slope_management_destroy(sm_plan)
    

#Need to ensure that we've initialized the numpy c-API.
np.import_array()

def slope_management_fast(np.ndarray[np.int_t, ndim=2] ap, np.ndarray[np.double_t, ndim=2] sy, np.ndarray[np.double_t, ndim=2] sx):
    """Slope management done in a single function call, including memory allocation, etc.
    
    This function acts in-place on sx and sy.
    """
    cdef int nx, ny
    ny = ap.shape[0]
    nx = ap.shape[1]
    slope_management(ny, nx, <int *>np.PyArray_DATA(ap.astype(np.int32)), <double *>np.PyArray_DATA(sy), <double *>np.PyArray_DATA(sx))
    return

cdef class SlopeManager:
    """A class to manage memory allocated for slope management on a specific aperture."""
    
    cdef sm_plan _plan
    cdef tuple _shape
    cdef np.ndarray _ap
    
    def __cinit__(self, ap):
        """C initialization step for this function."""
        shape = np.asarray(ap).shape
        self._shape = shape
        self._ap = ap
        
        # Create the plan, which internally allocates the necessary memory.
        self._plan = slope_management_plan(shape[0], shape[1], <int *>np.PyArray_DATA(ap.astype(np.int32)))
        
    @property
    def aperture(self):
        return self._ap
        
    def __dealloc__(self):
        slope_management_destroy(self._plan)
        
    def __call__(self, np.ndarray[np.double_t, ndim=2] sy, np.ndarray[np.double_t, ndim=2] sx):
        """Perform the actual slope management.
        
        The slope management happens in-place in the array.
        
        :parameter np.ndarray sy: The y slopes.
        :parameter np.ndarray sx: The x slopes.
        
        This function does not return anything because all work happens in place.
        """
        slope_management_execute(self._plan, <double *>np.PyArray_DATA(sy), <double *>np.PyArray_DATA(sx))
        return None
        