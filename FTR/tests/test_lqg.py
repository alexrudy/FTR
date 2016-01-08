# -*- coding: utf-8 -*-
"""
Tests for the LQG module.
"""

from __future__ import absolute_import
import pytest
import numpy as np

from .test_filter_base import FilterTestBase
from ..lqg import LQGFilter

class TestLQGFilter(FilterTestBase):
    """Tests for the LQG filter."""
    
    cls = LQGFilter
    
    repr = "<LQGFilter nlayers=... shape={shape}>"
    
    @pytest.fixture(params=[1,3,8])
    def nlayers(self, request):
        """Return the number of layers."""
        return request.param
    
    @pytest.fixture
    def filter(self, nlayers, shape):
        """Generate a random filter for testing."""
        full_shape = (nlayers,) + shape
        gains = np.random.randn(*full_shape) + 1j * np.random.randn(*full_shape)
        alphas = np.random.randn(*full_shape) + 1j * np.random.randn(*full_shape)
        hp_coeffs = np.random.randn(*shape) + 1j * np.random.randn(*shape)
        return self.cls(gains, alphas, hp_coeffs)
        
    