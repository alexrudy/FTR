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

def do_slope_management(np.ndarray[np.int_t, ndim=2] ap, np.ndarray[np.double_t, ndim=2] sy, np.ndarray[np.double_t, ndim=2] sx):
    """Slope management done in a single function call, including memory allocation, etc.
    
    This function acts in-place on sx and sy.
    """
    
    cdef int nx, ny
    ny = ap.shape[0]
    nx = ap.shape[1]
    _ap = ap.astype(np.int32)
    slope_management(ny, nx, <int *>np.PyArray_DATA(_ap), <double *>np.PyArray_DATA(sy), <double *>np.PyArray_DATA(sx))
    return
