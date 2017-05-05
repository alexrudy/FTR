# -*- coding: utf-8 -*-
"""
Base class for reconstructors, defines the reconstructor python API
"""

import abc
import six

@six.add_metaclass(abc.ABCMeta)
class Reconstructor(object):
    """A basic reconstructor."""
    
    @abc.abstractmethod
    def reconstruct(self, sx, sy):
        """Given the x and y slopes, reconstruct the phase.
        
        :param sx: x slopes.
        :param sy: y slopes.
        """
        pass
    
    @abc.abstractproperty
    def shape(self):
        """Shape of the reconstructed grid."""
        pass
        
    def __call__(self, xs, ys):
        """Reconstruct the phase.
        
        :param xs: The x slopes.
        :param ys: The y slopes.
        
        """
        return self.reconstruct(xs, ys)

@six.add_metaclass(abc.ABCMeta)
class Filter(object):
    """Base class for things which further filter phase, in the fourier domain."""
    
    def __call__(self, est_ft):
        """Callable to directly apply the filter."""
        return self.apply_filter(est_ft)
    
    @abc.abstractmethod
    def apply_filter(self, est_ft):
        """Apply the filter."""
        return est_ft