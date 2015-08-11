# -*- coding: utf-8 -*-
"""
A matrix-based reconstructor.
"""
from __future__ import (absolute_import, unicode_literals, division,
                        print_function)

# Python Imports
import six

# Scientific Python Imports
import numpy as np

# Local imports
from .base import Reconstructor
from .util import complexmp, ignoredivide

class MatrixReconstructor(Reconstructor):
    """A basic matrix reconstructor."""
    
    _shape = None
    """Internal storage for the shape of the reconstructor."""
    
    _matrix = None
    
    def __init__(self, matrix, shape=None):
        super(MatrixReconstructor, self).__init__()
        if shape is None:
            shape = np.asarray(matrix).shape[0:1]
        self._shape = tuple(shape)
        self.matrix = matrix
        
    @property
    def shape(self):
        """Shape of the phase input to this reconstructor."""
        return self._shape
        
    @property
    def matrix(self):
        """Ensure that the matrix works"""
        return self._matrix
        
    @matrix.setter
    def matrix(self, value):
        """Set the matrix."""
        matrix = np.asmatrix(value)
        # Check shape!
        if matrix.shape[0] != np.prod(self.shape):
            raise ValueError("Matrix output must be compatible with reconstructor shape. {!r} != {!r}".format(matrix.shape, self.shape))
        self._matrix = matrix
        
    def reconstruct(self, sx, sy):
        """Handle the matrix reconstruction."""
        s = np.concatenate((sx.flatten(), sy.flatten()))
        estimate = self.matrix * np.matrix([s]).T
        return np.asarray(estimate).reshape(self.shape)
        
