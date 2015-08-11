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
